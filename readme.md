# Orthogonal 12-mer barcodes for DropSynth beads

This repository contains barcode sets, source sequences, intermediate results, and historical analysis scripts for designing 12-nucleotide barcodes for DropSynth barcoded beads. Candidate sequences were derived from the Elledge 240,000-member 25-mer barcode library, screened for sequence composition, melting temperature, restriction sites, and self-dimerization, and pruned using Levenshtein distance.

The main bead barcode set contains **12,288 unique 12-mers**, selected as the first 12,288 entries of a 15,500-barcode set. A 6,144-barcode prefix subset is also included. The full set has no pair of sequences with Levenshtein distance below 3 in the stored orientation; both subsets inherit this property. This sequence-distance criterion does not by itself establish experimental hybridization specificity.

## Use the existing barcode sets

No software installation is needed to use the supplied FASTA files.

| File | Records | Purpose |
| --- | ---: | --- |
| [cp12mer\_15.5k_12288subset.fasta](cp12mer_15.5k_12288subset.fasta) | 12,288 | Main DropSynth bead barcode set |
| [cp12mer\_15.5k_6144subset.fasta](cp12mer_15.5k_6144subset.fasta) | 6,144 | Smaller bead barcode set |
| [cp12mer_15.5k.fasta](cp12mer_15.5k.fasta) | 15,500 | Full selected set |

Each file contains unique, unflanked 12-base DNA sequences with 50% GC content. Subset files preserve the sequence order and identifiers of the full set (`cp12mer_15.5k_1`, etc.). These are barcode sequences, not complete bead oligonucleotides.

## Repository contents

The original file locations are retained because the scripts use relative paths. Run Python and shell scripts from the repository root; open the R scripts in RStudio as described below.

| Files or directory | Role |
| --- | --- |
| `seq_files/` | Original 240,000 25-mers; SKPP15/SKPP20 forward and reverse sequences; older 384- and 1,536-barcode sets; corresponding BLAT `.2bit` files |
| `BCs_psl/` | Saved BLAT alignments (`.psl`) and readable alignments (`*_pretty.txt`); the historical `skpp20FWD_pretty.out` duplicates the corresponding `.txt` file |
| `BLAT_bcs.sh` | Aligns existing barcode sets against the original 25-mer library |
| `Generate_barcode_csv.py` | Parses readable BLAT alignments into a parent/subset mapping |
| `bc25mers_and_subsets_barcodes.csv` | Saved mapping, with parent names/sequences and matches to the six existing sets |
| `newBCs.R` | Selects parent 25-mers without SKPP15/SKPP20 matches |
| `unused_elledge_25mers.fasta` | 233,997 retained parent 25-mers |
| `1_12mer_barcodedesign_v5_noLev.py` | Generates and screens 12-mer candidates |
| `Primerselectiontools_py3.py` | Local melting-temperature and primer-dimer scoring functions |
| `restriction_minimal_Fall2021*.fasta`, `BspQI.fasta` | Restriction motifs used during candidate screening |
| `filt_prim_12nt_Tm_38_44_GC_45_55_SD_4.fasta` | Screened candidates: 189,935 records representing 184,146 unique sequences |
| `filt_lev.R` | Deduplicates candidates and prunes sequences by Levenshtein distance |
| `BCs202122.RData` | Historical saved R workspace; not explicitly loaded by the supplied scripts |
| `cp12mer_15.5k*.fasta` | Final barcode sets listed above |

## Software requirements for regeneration

These are historical scripts, rather than a packaged, version-locked pipeline. Exact original dependency versions were not recorded, and end-to-end regeneration has not been validated in a fresh environment.

- **UCSC tools:** `blat`, `faToTwoBit`, and `pslPretty`, available on `PATH`. See the [BLAT documentation](https://genome.ucsc.edu/goldenPath/help/blatSpec.html).
- **Python 3:** Biopython, pandas, matplotlib (provides `pylab`), and the `cLev` module. The latter is imported even though Python-side Levenshtein filtering is disabled by default; its implementation is not included here. The original README described it as “Cython cLev.”
- **R and RStudio:** `dplyr`, `tidyr`, `reshape2`, `magrittr`, `purrr`, `stringdist`, `Biostrings`, and `rstudioapi`. `Biostrings` is a Bioconductor package. Both R scripts use the active RStudio document to select their working directory.

Compatibility needs review before regeneration: the Python scripts use the legacy `Bio.SeqUtils.GC` API and positional arguments to pandas `str.split`. Use a compatible environment or adapt those calls before running. No tested environment lockfile is supplied.

## Historical workflow

Regeneration writes to the same filenames as the supplied results. Use a separate checkout or working copy to retain the archived outputs for comparison. Resolve the limitations below before treating a new run as a reproduction.

### 1. Map existing barcode sets to the parent library

```sh
bash BLAT_bcs.sh
python3 Generate_barcode_csv.py
```

The shell script aligns SKPP20 forward/reverse, SKPP15 forward/reverse, and the older 384- and 1,536-barcode sets against `seq_files/bc25mer.240k.fasta`. It writes alignments to `BCs_psl/`. The Python script parses the readable alignments and writes `bc25mers_and_subsets_barcodes.csv`.

The parser assumes two-line FASTA records and four-line readable alignment blocks. The saved mapping contains 240,005 rows; it is not a one-row-per-parent lookup, because alignment joins can create multiple rows per parent. Its repeated `loc_x`/`loc_y` columns are not emitted by the current parser, so the saved CSV and supplied parser differ in schema.

### 2. Select unused parent 25-mers

Open `newBCs.R` in RStudio and run it with that script active. It reads the mapping CSV and writes `unused_elledge_25mers.fasta`.

**The implemented selection excludes SKPP15 and SKPP20 matches only.** It does not filter the `bc384` or `bc1536` columns. In the saved mapping, 428 rows retained by this selection have a match to one or both older bead sets. Thus, “unused” in the filename should not be interpreted as exclusion of all six existing sets.

### 3. Screen 12-mer candidates

```sh
python3 1_12mer_barcodedesign_v5_noLev.py
```

The script scans parent 25-mers and keeps the first surviving candidate per parent after these filters:

| Criterion | Implemented setting |
| --- | --- |
| Barcode length | 12 nt |
| GC content | Strictly greater than 45% and less than 55%; exactly 50% for a 12-mer |
| Predicted melting temperature | 38–44 °C, inclusive, using `oligoTm` in `Primerselectiontools_py3.py` |
| Self-dimer score | At most 4, using the local `primerdimers` scoring function |
| Restriction motifs | BtsI, BspQI, BsmAI, BsrDI, BstNBI, BsaI, BbsI, BsmBI, BtgZI, Esp3I, and PaqCI, with reverse complements |
| Flanking context | Additional restriction checks on `GCTCTTCG` + barcode + `CGAAGAGC` |
| Python Levenshtein filtering | Disabled (`Lev_dist_check = False`); performed in R instead |

The BspQI-specific check permits no more than one occurrence of each BspQI motif in the flanked sequence. `BspQI.fasta` contains the two motifs copied from `restriction_minimal_Fall2021.fasta`; this missing input was restored during repository preparation.

The candidate output is `filt_prim_12nt_Tm_38_44_GC_45_55_SD_4.fasta`, containing **189,935 records / 184,146 unique sequences** in the archived file. The script also produces melting-temperature histograms.

The historical extraction loop uses `range(25 - 12)`: it considers 13 windows per parent and omits the final possible 12-mer window. Changing this behavior would change the candidate pool.

### 4. Prune candidates by Levenshtein distance

Open `filt_lev.R` in RStudio and run it with that script active. It deduplicates the candidate sequences, prunes pairs with Levenshtein distance below 3 within chunks, and then prunes merged groups. It writes `cp12mer_15.5k.fasta`; the supplied result contains 15,500 sequences.

**The supplied script processes only the first 160,000 unique candidates**, in eight chunks of 20,000. It does not process the remaining 24,146 unique candidates in the archived input. The greedy selection depends on candidate order and chunk boundaries.

This step is memory intensive: one 20,000-sequence chunk creates 199,990,000 pairs before filtering, and merged groups also require pairwise comparisons. Chunking does not make this a low-memory operation.

### 5. Select the bead barcode subsets

The saved 12,288- and 6,144-barcode sets are exact prefixes of the full set. Their export is not included in `filt_lev.R`. Once `cp12mer_15.5k.fasta` exists, they can be exported in R using:

```r
library(Biostrings)
barcodes <- readDNAStringSet("cp12mer_15.5k.fasta")
stopifnot(length(barcodes) >= 12288L)
writeXStringSet(barcodes[seq_len(12288)], "cp12mer_15.5k_12288subset.fasta")
writeXStringSet(barcodes[seq_len(6144)], "cp12mer_15.5k_6144subset.fasta")
```

## Verification of the supplied sets

Repository checks during publication preparation confirmed:

- The full set contains 15,500 unique 12-mers, all with six G/C bases, and all occur in the candidate FASTA.
- The 12,288- and 6,144-member sets exactly match the corresponding prefixes of the full set, including identifiers.
- No two sequences in the full set have Levenshtein distance below 3. This was checked independently by enumerating one- and two-substitution neighbors and shared single-deletion signatures (covering insertion/deletion pairs for equal-length sequences).

These checks apply to the supplied sequence files in their stored orientation. They do not constitute an end-to-end pipeline rerun or a reverse-complement distance check.

## Sources and license

The parent barcode library comes from the Elledge laboratory's *Design of 240,000 orthogonal 25mer DNA barcode probes*; see the [Elledge barcode resource](https://elledge.hms.harvard.edu/?page_id=638). The alignment-parsing code was attributed in the original README to Nora's work in [PlesaLab/Barcodes\_from\_Elledge](https://github.com/PlesaLab/Barcodes_from_Elledge).

Repository code is distributed under the [MIT license](LICENSE). The source links above document the provenance of the external barcode library and alignment parser.
