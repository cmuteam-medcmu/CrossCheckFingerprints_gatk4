# CrossCheckFingerprints_gatk4

A Python wrapper around **GATK4 / Picard `CrosscheckFingerprints`** for checking sample identity across many
sequencing files at once.

Given a list of BAM/CRAM/VCF files, the tool fingerprints each file at a set of common SNP sites, compares every
file against every other file, and reports:

1. **`cross_all_sample.csv`** — the full pairwise comparison table (LOD scores and MATCH / MISMATCH results).
2. **`inferred_group.csv`** — files clustered into inferred donors/individuals, so you can immediately see which
   files came from the same person.

Typical uses:

- Detecting **sample swaps or mislabeling** in a sequencing batch.
- Confirming that **tumor/normal pairs** really come from the same patient.
- Checking that multiple runs, lanes, or assays of the "same" sample are in fact the same individual.
- QC of large cohorts, where running `CrosscheckFingerprints` by hand and reading the raw metrics file is
  impractical.

---

## How it works

```
sample_full_path.txt        haplotype map (.txt)
        │                          │
        └──────────┬───────────────┘
                   ▼
        gatk CrosscheckFingerprints
         (--CROSSCHECK_BY FILE)
                   │
                   ▼
        raw crosscheck metrics
                   │
        ┌──────────┴───────────┐
        ▼                      ▼
 cross_all_sample.csv    inferred_group.csv
 (all pairwise LODs)     (files grouped by donor)
```

The LOD score is the log-odds that two files come from the same individual. Positive LOD → same individual,
negative LOD → different individuals. Grouping is done by linking together all file pairs whose LOD exceeds the
match threshold.

---

## Requirements

| Component | Notes |
|---|---|
| **GATK 4** | Tested with `gatk/4.6.2.0`. Must be on `$PATH` when running with `--no-docker`. |
| **Java 17** | Required by GATK 4.6. Usually loaded together with the GATK module. |
| **Python 3** | Python 3.8 or newer. |
| **pysam** | `pip install pysam` |
| **Docker** | Only needed if you *don't* pass `--no-docker` (runs GATK from a container image instead of a local install). |

On an HPC system with environment modules:

```bash
ml purge
ml gatk/4.6.2.0
pip install --user pysam
```

### Haplotype map

You also need a **haplotype database** matching your reference build, e.g. for GRCh38:

```
/path/to/hg38/v0/Homo_sapiens_assembly38.haplotype_database.txt
```

This ships with the GATK resource bundle. Pre-built maps for other builds and for RNA-seq / ATAC-seq style data
are available at <https://github.com/naumanjaved/fingerprint_maps>.

> The contig names and order in the haplotype map **must match the reference used to align your BAMs**
> (e.g. `chr1` vs `1`), otherwise GATK will error out.

---

## Installation

```bash
git clone https://github.com/cmuteam-medcmu/CrossCheckFingerprints_gatk4.git
cd CrossCheckFingerprints_gatk4
pip install pysam
```

---

## Usage

### 1. Prepare the input list

Create `workdir/sample_full_path.txt` with one absolute file path per line:

```
/data/project/align/S001.sorted.bam
/data/project/align/S002.sorted.bam
/data/project/align/S003.sorted.bam
```

Each BAM should be coordinate-sorted and indexed (`.bai` next to the file).

### 2. Run

Edit and run the provided wrapper script:

```bash
bash run_script.sh
```

`run_script.sh` contains:

```bash
#!/bin/bash
# ml purge
# ml gatk/4.6.2.0
# python3 requires: pip install pysam

CROSS_CHECK='src/crosscheck_fingerprints.py'
WORKDIR='workdir'
PANEL='/path/to/hg38/v0/Homo_sapiens_assembly38.haplotype_database.txt'

python3 ${CROSS_CHECK} \
    --input-list    ${WORKDIR}/sample_full_path.txt \
    --haplotype-map ${PANEL} \
    --output        cross_all_sample.csv \
    --group-output  inferred_group.csv \
    --workdir       ${WORKDIR} \
    --crosscheck-by FILE \
    --id-suffix     "sorted.bam" \
    --no-docker \
    --num-threads   2
```

Before the first run, set `PANEL` to your local haplotype map.

---

## Options

| Option | Description |
|---|---|
| `--input-list` | Text file with one input file path per line (BAM / CRAM / VCF). |
| `--haplotype-map` | Haplotype database for the matching reference build. |
| `--output` | Output CSV with all pairwise comparisons. Default: `cross_all_sample.csv`. |
| `--group-output` | Output CSV with files clustered into inferred individuals. |
| `--workdir` | Directory for intermediate files and logs. |
| `--crosscheck-by` | Comparison unit: `FILE`, `SAMPLE`, `LIBRARY`, or `READGROUP`. Use `FILE` when read-group metadata is unreliable. |
| `--id-suffix` | Suffix stripped from filenames to derive the sample ID (e.g. `sorted.bam` turns `S001.sorted.bam` into `S001`). |
| `--no-docker` | Use a locally installed `gatk` instead of running it inside Docker. |
| `--num-threads` | Threads passed to GATK (`NUM_THREADS`). |

Run `python3 src/crosscheck_fingerprints.py --help` for the authoritative list.

---

## Output

### `cross_all_sample.csv`

One row per file pair:

| Column | Meaning |
|---|---|
| `LEFT_FILE` / `RIGHT_FILE` | The two files being compared. |
| `RESULT` | `MATCH`, `MISMATCH`, or `INCONCLUSIVE`. |
| `LOD_SCORE` | Log-odds that the two files are from the same individual. |
| `LOD_SCORE_TUMOR_NORMAL` | Tumor-aware LOD, robust to loss of heterozygosity. |

### `inferred_group.csv`

Each input file with an assigned **INFERRED_GROUP** ID. Files sharing a **INFERRED_GROUP** ID are inferred to come from the same
individual — compare this against your sample sheet to spot swaps.

---

## Interpreting the results

- **LOD > 5** — strong evidence the two files are the same individual.
- **LOD < -5** — strong evidence they are different individuals.
- **-5 to 5** — inconclusive; usually too few reads at fingerprint sites.

Low coverage, targeted panels, and RNA-seq all reduce the number of usable fingerprint sites and push results
toward inconclusive. If most comparisons are inconclusive, check that (a) the haplotype map matches your
reference, and (b) your data actually covers the fingerprint sites.
