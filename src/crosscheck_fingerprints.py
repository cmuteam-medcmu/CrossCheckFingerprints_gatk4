#!/usr/bin/env python3
"""
crosscheck_fingerprints.py

Wraps GATK4's CrosscheckFingerprints (Picard) to verify sample identity by
comparing genotypes at a curated fingerprinting SNP panel. Accepts BAM, CRAM,
or VCF inputs.

Uses pysam to validate inputs (index present, sample names in header) before
invoking GATK, then converts the Picard metrics report into a condensed CSV
with clean sample IDs (derived from filenames) and the key result columns
(RESULT, DATA_TYPE, LOD scores, lane/barcode/library). Pass --raw-output to
get GATK's full, unmodified metrics table instead.

Docs: https://gatk.broadinstitute.org/hc/en-us/articles/360037594711-CrosscheckFingerprints-Picard-

Requirements:
    pip install pysam
    GATK4 available either as a local `gatk` executable on PATH, or via Docker
    (default image: broadinstitute/gatk:4.6.1.0)

Examples:
    # Compare a BAM against a second BAM (e.g. expected/prior sample)
    python3 crosscheck_fingerprints.py \\
        --input recalibrated.bam \\
        --second-input expected_sample.bam \\
        --haplotype-map hg38_chr.map \\
        --output crosscheck_results.csv

    # Cross-check all samples within a single VCF against each other
    python3 crosscheck_fingerprints.py \\
        --input output.vcf \\
        --haplotype-map hg38_chr.map \\
        --output crosscheck_results.csv

    # Cross-check many samples using a text file of paths (one per line);
    # file lists can be combined with each other and with -1/-2 directly
    python3 crosscheck_fingerprints.py \\
        --input-list cohort_a.txt \\
        --second-input-list cohort_b.txt \\
        --haplotype-map hg38_chr.map \\
        --output crosscheck_results.csv

    # Use a local gatk install instead of docker
    python3 crosscheck_fingerprints.py \\
        --input recalibrated.bam \\
        --haplotype-map hg38_chr.map \\
        --output crosscheck_results.csv \\
        --no-docker
"""

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import pysam
except ImportError:
    pysam = None


CONTAINER_WORKDIR = "/workdir"

# Columns kept in the summarized output, in order: (output_name, source_name)
SUMMARY_COLUMNS = [
    ("RESULT", "RESULT"),
    ("DATA_TYPE", "DATA_TYPE"),
    ("LOD", "LOD_SCORE"),
    ("LOD_TUMOR_NORMAL", "LOD_SCORE_TUMOR_NORMAL"),
    ("LOD_NORMAL_TUMOR", "LOD_SCORE_NORMAL_TUMOR"),
    ("LEFT_LANE", "LEFT_LANE"),
    ("RIGHT_LANE", "RIGHT_LANE"),
    ("LEFT_BARCODE", "LEFT_RUN_BARCODE"),
    ("RIGHT_BARCODE", "RIGHT_RUN_BARCODE"),
    ("LEFT_LIBRARY", "LEFT_LIBRARY"),
    ("RIGHT_LIBRARY", "RIGHT_LIBRARY"),
]


def parse_args():
    p = argparse.ArgumentParser(
        description="Run GATK CrosscheckFingerprints on BAM/CRAM/VCF input(s) "
                    "and save the results as CSV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-1", "--input", action="append", dest="inputs",
                   help="Sample 1: input file (BAM/CRAM/VCF). Repeatable; combines "
                        "with --input-list.")
    p.add_argument("--input-list",
                   help="Text file listing --input paths, one per line. Blank lines "
                        "and lines starting with '#' are ignored.")
    p.add_argument("-2", "--second-input", action="append", dest="second_inputs",
                   help="Sample 2: file(s) to cross-check --input against. If omitted, "
                        "samples within --input are cross-checked against each other. "
                        "Repeatable; combines with --second-input-list.")
    p.add_argument("--second-input-list",
                   help="Text file listing --second-input paths, one per line.")
    p.add_argument("--haplotype-map",
                   default="/common/db/human_ref/hg38/v0/Homo_sapiens_assembly38.haplotype_database.txt",
                   help="Haplotype map file defining the fingerprint SNP panel.")
    p.add_argument("--reference",
                   help="Reference FASTA (required for CRAM, recommended for VCF).")
    p.add_argument("--output", required=True, help="Path to write the final CSV results.")
    p.add_argument("--group-output",
                   help="Also write a CSV report inferring which individual each "
                        "sample belongs to. Columns: ID, INFERRED_GROUP, MIN_LOD, "
                        "MEAN_LOD, MAX_LOD, NOTE. One row per unique LEFT_ID sample; "
                        "LOD stats are computed across every raw comparison row for "
                        "that sample (no rows are dropped/deduped for this report). "
                        "Samples connected by a MATCH (directly or transitively) share "
                        "an INFERRED_GROUP code (00_A, 00_B, ... 00_Z, 01_A, ...; sorts "
                        "in generation order) and have a blank NOTE; MIN/MEAN/MAX_LOD "
                        "are computed only from that sample's comparisons whose "
                        "RIGHT_ID also lands in the same inferred group. Samples with "
                        "no MATCH at all are inferred from their top-2 positive-LOD "
                        "comparison partners: if both partners share a group, this "
                        "sample is folded into it and noted 'LOW_LOD|POSITIVE_INFERRED'; "
                        "otherwise it gets its own group and is noted 'ALL INCLUSIVE'. "
                        "Any sample left as the sole member of its group in this report "
                        "is also forced to 'ALL INCLUSIVE' with '-' in the LOD columns.")
    p.add_argument("--crosscheck-mode", default="CHECK_ALL_OTHERS",
                   choices=["CHECK_SAME_SAMPLE", "CHECK_ALL_OTHERS"],
                   help="GATK CROSSCHECK_MODE. CHECK_ALL_OTHERS (default) compares every "
                        "sample regardless of name -- use this unless sample names are "
                        "guaranteed to match across files. CHECK_SAME_SAMPLE only compares "
                        "identically-named samples (requires --crosscheck-by SAMPLE) and "
                        "reports unmatched ones instead of comparing them.")
    p.add_argument("--crosscheck-by", default="SAMPLE",
                   choices=["SAMPLE", "FILE", "READGROUP", "LIBRARY"],
                   help="Grouping level for comparison (default: SAMPLE).")
    p.add_argument("--lod-threshold", type=float, default=-5.0,
                   help="LOD threshold for match/mismatch calls (default: -5.0).")
    p.add_argument("--exit-code-when-mismatch", type=int, default=1,
                   help="Exit code GATK returns on mismatch (default: 1). The CSV is "
                        "still written.")
    p.add_argument("--num-threads", type=int, default=4,
                   help="Threads for GATK to use (default: 4).")
    p.add_argument("--workdir",
                   help="Directory to mount/use for inputs (default: common parent of "
                        "all input files).")
    p.add_argument("--gatk-image", default="broadinstitute/gatk:4.6.1.0",
                   help="Docker image for GATK (default: broadinstitute/gatk:4.6.1.0).")
    p.add_argument("--no-docker", action="store_true",
                   help="Use a local `gatk` executable on PATH instead of Docker.")
    p.add_argument("--skip-validation", action="store_true",
                   help="Skip pysam-based pre-flight validation of inputs.")
    p.add_argument("--keep-raw-metrics", action="store_true",
                   help="Also keep the raw Picard metrics file next to the CSV.")
    p.add_argument("--id-suffix", default="_recal.bam",
                   help="Suffix stripped from LEFT_FILE/RIGHT_FILE basenames to derive "
                        "LEFT_ID/RIGHT_ID (default: '_recal.bam'), e.g. "
                        "'OS059_Dx_1_recal.bam' -> 'OS059_Dx_1'.")
    p.add_argument("--raw-output", action="store_true",
                   help="Write GATK's full metrics table instead of the condensed summary.")
    p.add_argument("--keep-self-comparisons", action="store_true",
                   help="Keep rows where LEFT_ID == RIGHT_ID (dropped by default).")
    p.add_argument("--keep-duplicate-pairs", action="store_true",
                   help="Keep both A-vs-B and B-vs-A (only one is kept by default).")
    p.add_argument("--keep-full-result", action="store_true",
                   help="Keep GATK's full RESULT value. By default the EXPECTED_/"
                        "UNEXPECTED_ prefix is dropped, leaving MATCH or MISMATCH.")

    args = p.parse_args()

    # Merge file-list paths with paths given directly via -1/-2.
    args.inputs = (args.inputs or []) + read_file_list(args.input_list, "--input-list")
    args.second_inputs = (args.second_inputs or []) + read_file_list(
        args.second_input_list, "--second-input-list")

    if not args.inputs:
        p.error("no input files given: provide -1/--input and/or --input-list")

    return args


def read_file_list(list_path, flag_name):
    """Read a text file of paths, one per line, ignoring blanks and '#' comments."""
    if not list_path:
        return []
    if not Path(list_path).exists():
        raise FileNotFoundError(f"{flag_name} file not found: {list_path}")

    paths = [line.strip() for line in Path(list_path).read_text().splitlines()]
    paths = [line for line in paths if line and not line.startswith("#")]
    if not paths:
        print(f"!!! WARNING: {flag_name} file '{list_path}' contained no paths.",
              file=sys.stderr)
    return paths


def validate_inputs(paths, reference=None):
    """Check that alignment files are indexed and readable, and log the sample
    names found in each header."""
    if pysam is None:
        print("!!! pysam is not installed - skipping input validation "
              "(pip install pysam to enable).", file=sys.stderr)
        return

    for path in paths:
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"Input file not found: {path}")
        ext = p.suffix.lower()

        try:
            if ext in (".bam", ".cram"):
                kwargs = {"reference_filename": reference} if ext == ".cram" else {}
                mode = "rc" if ext == ".cram" else "rb"
                with pysam.AlignmentFile(str(p), mode, **kwargs) as af:
                    samples = sorted({rg.get("SM", "UNKNOWN")
                                      for rg in af.header.get("RG", [])})
                    indexed = af.has_index()
                if not indexed:
                    print(f"    [validate] WARNING: {p.name} has no index. "
                          f"CrosscheckFingerprints requires a .bai/.crai.", file=sys.stderr)
            elif ext == ".vcf" or p.name.endswith(".vcf.gz"):
                with pysam.VariantFile(str(p)) as vf:
                    samples = list(vf.header.samples)
            else:
                print(f"    [validate] {p.name}: unrecognized extension, skipping deep check.")
                continue
            print(f"    [validate] {p.name}: samples in header -> "
                  f"{', '.join(samples) or '(none found)'}")
        except Exception as e:
            print(f"    [validate] WARNING: could not fully validate {p.name}: {e}",
                  file=sys.stderr)


def build_gatk_command(args, workdir, metrics_path):
    """Build the CrosscheckFingerprints command, either as a local `gatk`
    invocation or wrapped in `docker run`."""

    def arg_path(local_path):
        p = Path(local_path).resolve()
        if args.no_docker:
            return str(p)
        try:
            return str(Path(CONTAINER_WORKDIR) / p.relative_to(workdir))
        except ValueError:
            raise SystemExit(
                f"!!! {p} is outside the mounted workdir ({workdir}); docker cannot see "
                f"it. Pass a --workdir that contains every input, or use --no-docker.")

    cmd = ["gatk", "CrosscheckFingerprints"]
    for inp in args.inputs:
        cmd += ["--INPUT", arg_path(inp)]
    for sinp in args.second_inputs:
        cmd += ["--SECOND_INPUT", arg_path(sinp)]
    cmd += ["--HAPLOTYPE_MAP", arg_path(args.haplotype_map)]
    if args.reference:
        cmd += ["-R", arg_path(args.reference)]
    cmd += [
        "--OUTPUT", arg_path(metrics_path),
        "--CROSSCHECK_BY", args.crosscheck_by,
        "--CROSSCHECK_MODE", args.crosscheck_mode,
        "--LOD_THRESHOLD", str(args.lod_threshold),
        "--EXIT_CODE_WHEN_MISMATCH", str(args.exit_code_when_mismatch),
        "--NUM_THREADS", str(args.num_threads),
    ]

    if args.no_docker:
        return cmd
    return ["docker", "run", "--rm",
            "--volume", f"{workdir}:{CONTAINER_WORKDIR}",
            "--workdir", CONTAINER_WORKDIR,
            args.gatk_image] + cmd


def run_gatk(cmd):
    print(">>> Running: " + " ".join(cmd))
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError as e:
        print(f"!!! Could not execute '{cmd[0]}': {e}", file=sys.stderr)
        print("!!! If using --no-docker, make sure `gatk` is on PATH "
              "(e.g. `module load gatk4`) in this job's environment.", file=sys.stderr)
        raise
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    return result


def parse_picard_metrics(metrics_file, stdout="", stderr=""):
    """Parse a Picard-style metrics file (comment header, then '## METRICS CLASS'
    section with a column header + rows) into (columns, rows)."""
    metrics_file = Path(metrics_file)
    if not metrics_file.exists():
        def tail(text):
            return "\n".join(text.strip().splitlines()[-25:])
        raise FileNotFoundError(
            f"Expected metrics file was not produced: {metrics_file}\n"
            "The GATK command failed before writing output. Full GATK output is above; "
            "here is what was captured:\n"
            f"--- GATK stderr (tail) ---\n{tail(stderr)}\n"
            f"--- GATK stdout (tail) ---\n{tail(stdout)}")

    lines = metrics_file.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith("## METRICS CLASS"):
            break
    else:
        raise ValueError(
            f"Could not find '## METRICS CLASS' section in {metrics_file}. "
            f"The GATK run may have failed before producing results.")

    columns = lines[i + 1].split("\t")
    rows = []
    for line in lines[i + 2:]:
        if not line.strip():
            break
        rows.append(line.split("\t"))
    return columns, rows


def derive_id(file_uri, suffix):
    """'file:///path/to/OS059_Dx_1_recal.bam' -> 'OS059_Dx_1'. Falls back to the
    file stem if the suffix isn't present."""
    if not file_uri:
        return file_uri
    name = Path(re.sub(r"^file://", "", file_uri)).name
    return name[:-len(suffix)] if suffix and name.endswith(suffix) else Path(name).stem


def summarize_rows(columns, rows, args):
    """Condense the raw columns/rows to LEFT_ID, RIGHT_ID + SUMMARY_COLUMNS,
    dropping self-comparisons and reversed duplicate pairs unless asked to keep
    them, and trimming RESULT to MATCH/MISMATCH."""
    required = {"LEFT_FILE", "RIGHT_FILE"} | {src for _, src in SUMMARY_COLUMNS}
    missing = required - set(columns)
    if missing:
        raise ValueError(
            f"Cannot summarize: metrics output is missing expected column(s) "
            f"{sorted(missing)}. Re-run with --raw-output to get the full table.")

    idx = {name: columns.index(name) for name in required}
    out_columns = ["LEFT_ID", "RIGHT_ID"] + [name for name, _ in SUMMARY_COLUMNS]
    result_i = out_columns.index("RESULT")

    out_rows, seen = [], set()
    for row in rows:
        left = derive_id(row[idx["LEFT_FILE"]], args.id_suffix)
        right = derive_id(row[idx["RIGHT_FILE"]], args.id_suffix)

        if left == right and not args.keep_self_comparisons:
            continue
        if not args.keep_duplicate_pairs:
            pair = frozenset((left, right))
            if pair in seen:
                continue
            seen.add(pair)

        new_row = [left, right] + [row[idx[src]] for _, src in SUMMARY_COLUMNS]
        if not args.keep_full_result:
            new_row[result_i] = re.sub(r"^(UN)?EXPECTED_", "", new_row[result_i])
        out_rows.append(new_row)

    return out_columns, out_rows


def generate_group_codes(count):
    """Return `count` short group codes, in generation order, formatted so
    that sorting them as plain strings reproduces that same order (like a
    timestamp-style ID). The digit block comes first and is zero-padded to
    a fixed width, the letter comes second: 00_A, 00_B, ..., 00_Z, 01_A,
    01_B, ... -- so code #26 (01_A) sorts immediately after code #25 (00_Z)
    instead of two codes that are far apart in generation order (e.g. the
    old A_0/A_1 scheme) accidentally sorting next to each other just
    because they share a leading letter.

    The digit width is sized to `count` so it never needs more than 2
    digits' worth of headroom, keeping string-sort correct at any scale.
    """
    letters = [chr(c) for c in range(ord("A"), ord("Z") + 1)]
    max_digit = (max(count - 1, 0)) // len(letters)
    digit_width = max(2, len(str(max_digit)))
    return [
        f"{n // len(letters):0{digit_width}d}_{letters[n % len(letters)]}"
        for n in range(count)
    ]


def assign_group_codes(raw_columns, raw_rows, id_suffix):
    """Union-find every sample ID over MATCH edges only (a MATCH row unions
    its LEFT_ID and RIGHT_ID; anything else leaves them apart). Returns
    (find, code_by_root, all_ids) where code_by_root maps each connected
    component's root to a short, sortable code from generate_group_codes(),
    assigned in a deterministic order (sorted by the component's smallest
    member ID)."""
    lf_i, rf_i, res_i = (raw_columns.index(c) for c in ("LEFT_FILE", "RIGHT_FILE", "RESULT"))

    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]  # path compression
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    all_ids = set()
    for row in raw_rows:
        left = derive_id(row[lf_i], id_suffix)
        right = derive_id(row[rf_i], id_suffix)
        if not left or not right:
            continue
        find(left)
        find(right)
        all_ids.add(left)
        all_ids.add(right)
        result = re.sub(r"^(UN)?EXPECTED_", "", row[res_i])
        if result == "MATCH":
            union(left, right)

    components = {}
    for sample_id in all_ids:
        components.setdefault(find(sample_id), []).append(sample_id)

    ordered_roots = sorted(components, key=lambda r: min(components[r]))
    codes = generate_group_codes(len(ordered_roots))
    code_by_root = dict(zip(ordered_roots, codes))
    return find, code_by_root, all_ids


def infer_group_rows(raw_columns, raw_rows, id_suffix):
    """Build the per-sample group-inference report: one row per unique
    LEFT_ID, with columns ID, INFERRED_GROUP, MIN_LOD, MEAN_LOD, MAX_LOD,
    NOTE. Runs on the raw (unfiltered, undeduplicated) metrics rows so every
    comparison is available -- nothing is dropped the way summarize_rows()
    drops self-comparisons/duplicate pairs.

    Grouping rules:
      1. If a sample has at least one MATCH result (as LEFT_ID or RIGHT_ID),
         its group is the connected component of MATCH edges it belongs to
         (see assign_group_codes). NOTE is left blank.
      2. Otherwise the sample is "inconclusive" (no MATCH at all). Look only
         at its comparisons where it is LEFT_ID and LOD is positive, take
         the top 2 by LOD:
           - if those two partner samples share a group, fold this sample
             into that group and set NOTE to "LOW_LOD|POSITIVE_INFERRED".
           - otherwise (fewer than 2 positive-LOD partners, or the top two
             land in different groups), the sample keeps its own singleton
             group and NOTE is "ALL INCLUSIVE".

    Once a sample's group is settled, MIN_LOD/MEAN_LOD/MAX_LOD are computed
    only from that sample's (LEFT_ID) comparisons whose RIGHT_ID partner
    also belongs to that same inferred group -- not from every comparison.
    """
    lf_i, rf_i, res_i = (raw_columns.index(c) for c in ("LEFT_FILE", "RIGHT_FILE", "RESULT"))
    lod_col = "LOD" if "LOD" in raw_columns else "LOD_SCORE"
    lod_i = raw_columns.index(lod_col) if lod_col in raw_columns else None

    find, code_by_root, _ = assign_group_codes(raw_columns, raw_rows, id_suffix)

    comparisons_by_id = {}   # left_id -> [(right_id, lod_or_None), ...]  (left_id is LEFT_ID)
    has_match = {}           # sample_id -> True if involved in >=1 MATCH row

    for row in raw_rows:
        left = derive_id(row[lf_i], id_suffix)
        right = derive_id(row[rf_i], id_suffix)
        if not left or not right:
            continue
        result = re.sub(r"^(UN)?EXPECTED_", "", row[res_i])

        has_match.setdefault(left, False)
        has_match.setdefault(right, False)
        if result == "MATCH":
            has_match[left] = True
            has_match[right] = True

        lod = None
        if lod_i is not None:
            try:
                lod = float(row[lod_i])
            except (ValueError, IndexError):
                lod = None

        comparisons_by_id.setdefault(left, [])
        comparisons_by_id[left].append((right, lod))

    out_rows = []
    for sample_id in sorted(comparisons_by_id):
        comparisons = comparisons_by_id[sample_id]

        if has_match.get(sample_id):
            group = code_by_root[find(sample_id)]
            note = ""
        else:
            positive = sorted(
                ((r, l) for r, l in comparisons if l is not None and l > 0),
                key=lambda pair: pair[1], reverse=True,
            )
            top2 = positive[:2]
            if len(top2) == 2 and code_by_root[find(top2[0][0])] == code_by_root[find(top2[1][0])]:
                group = code_by_root[find(top2[0][0])]
                note = "LOW_LOD|POSITIVE_INFERRED"
            else:
                group = code_by_root[find(sample_id)]
                note = "ALL INCLUSIVE"

        # MIN/MEAN/MAX only over comparisons whose RIGHT_ID lands in this
        # sample's inferred group (i.e. the partner's own group code matches).
        in_group_lods = [
            l for r, l in comparisons
            if l is not None and code_by_root.get(find(r)) == group
        ]
        min_lod = min(in_group_lods) if in_group_lods else None
        max_lod = max(in_group_lods) if in_group_lods else None
        mean_lod = (sum(in_group_lods) / len(in_group_lods)) if in_group_lods else None

        out_rows.append([
            sample_id,
            group,
            f"{min_lod:.4f}" if min_lod is not None else "",
            f"{mean_lod:.4f}" if mean_lod is not None else "",
            f"{max_lod:.4f}" if max_lod is not None else "",
            note,
        ])

    # If a sample ends up the only member of its inferred group in this
    # report, that's a group of one -- force NOTE to "ALL INCLUSIVE" and the
    # LOD columns to "-" regardless of how the group was originally settled.
    group_counts = {}
    for row in out_rows:
        group_counts[row[1]] = group_counts.get(row[1], 0) + 1
    for row in out_rows:
        if group_counts[row[1]] == 1:
            row[2] = row[3] = row[4] = "-"
            row[5] = "ALL INCLUSIVE"

    columns = ["ID", "INFERRED_GROUP", "MIN_LOD", "MEAN_LOD", "MAX_LOD", "NOTE"]
    return columns, out_rows


def write_groups(columns, rows, path):
    """Write the per-sample group-inference report as CSV."""
    write_csv(columns, rows, path)


def write_csv(columns, rows, csv_file):
    with open(csv_file, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(columns)
        writer.writerows(rows)


def report_mismatches(columns, rows, lod_threshold):
    lod_col = "LOD" if "LOD" in columns else "LOD_SCORE"
    if lod_col not in columns:
        return
    lod_i = columns.index(lod_col)
    for row in rows:
        try:
            lod = float(row[lod_i])
        except (ValueError, IndexError):
            continue
        flag = "MISMATCH" if lod < lod_threshold else "match"
        print(f"    LOD={lod:.2f} -> {flag}", file=sys.stderr)


def main():
    args = parse_args()

    all_paths = args.inputs + args.second_inputs + [args.haplotype_map]
    if args.reference:
        all_paths.append(args.reference)

    print(f">>> {len(args.inputs)} input file(s)"
          + (f", {len(args.second_inputs)} second-input file(s)" if args.second_inputs else ""))

    if not args.skip_validation:
        print(">>> Validating inputs ...")
        validate_inputs(args.inputs + args.second_inputs, reference=args.reference)

    if args.workdir:
        workdir = Path(args.workdir).resolve()
    else:
        workdir = Path(os.path.commonpath([str(Path(p).resolve().parent) for p in all_paths]))

    out_csv = Path(args.output).resolve()
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    # The metrics file must live under workdir so it resolves the same way
    # inside the container as it does on the host.
    metrics_path = workdir / f".crosscheck_metrics_{os.getpid()}.txt"
    result = run_gatk(build_gatk_command(args, workdir, metrics_path))

    try:
        raw_columns, raw_rows = parse_picard_metrics(metrics_path, result.stdout, result.stderr)
        if args.raw_output:
            columns, rows = raw_columns, raw_rows
        else:
            columns, rows = summarize_rows(raw_columns, raw_rows, args)
        write_csv(columns, rows, out_csv)

        if args.group_output:
            group_columns, group_rows = infer_group_rows(raw_columns, raw_rows, args.id_suffix)
            group_path = Path(args.group_output).resolve()
            group_path.parent.mkdir(parents=True, exist_ok=True)
            write_groups(group_columns, group_rows, group_path)

            n_groups = len({row[1] for row in group_rows})
            n_inferred = sum(1 for row in group_rows if row[-1])
            print(f">>> {n_groups} inferred group(s) across {len(group_rows)} sample(s), "
                  f"written to: {group_path}")
            if n_inferred:
                print(f"!!! {n_inferred} sample(s) had no direct MATCH and were grouped "
                      f"by inference instead -- check the NOTE column in {group_path}.",
                      file=sys.stderr)
    finally:
        if metrics_path.exists():
            if args.keep_raw_metrics:
                kept = out_csv.with_suffix(".picard_metrics.txt")
                shutil.copy(metrics_path, kept)
                print(f">>> Raw Picard metrics kept at: {kept}")
            metrics_path.unlink()

    print(f">>> CSV results written to: {out_csv} ({len(rows)} row(s))")

    if result.returncode != 0:
        print(f"!!! GATK exited with status {result.returncode}. This likely indicates "
              f"a fingerprint MISMATCH (LOD below threshold) or an inconclusive/"
              f"expected-mismatch result.", file=sys.stderr)
        report_mismatches(columns, rows, args.lod_threshold)

    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
