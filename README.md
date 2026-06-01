# Chimera — Structural Variant Fusion Gene Pipeline

Chimera calls candidate gene fusions from structural variant (SV) data, classifies them using an ACMG/AMP-inspired scoring framework, and produces a self-contained interactive HTML report with embedded SVG diagrams.

---

## Features

- **Two input formats** — headerless BED-style TSV or MAVIS-style headered TSV with breakpoint orientation columns
- **Exact breakpoint–gene overlap** using a sorted interval index (no window expansion, no spurious hits)
- **Canonical fusion classification** using a strand × SV-type decision matrix covering deletions, duplications, inversions, translocations, and inverted translocations
- **Reading frame prediction** in superTranscript space using GENCODE CDS records, preferring MANE Select / Ensembl canonical transcripts
- **ACMG/AMP-inspired pathogenicity scoring** integrating ClinGen dosage sensitivity scores and COSMIC Cancer Gene Census
- **Interactive HTML report** — tabbed interface with cohort summary, fusion cards with inline SVG diagrams, raw data tables, visualisations (circos plot, rank plot, force-directed network), and an IGV/UCSC integration guide
- **IGV-ready outputs** — BED file of breakpoints colour-coded by pathogenicity, plus a pre-configured IGV session XML

---

## Requirements

- Python ≥ 3.9
- `pandas`

```bash
pip install pandas
```

---

## Installation

```bash
git clone https://github.com/your-org/chimera.git
cd chimera
pip install pandas
```

No build step is required — the pipeline is a single Python script.

---

## Quick Start

```bash
python chimera.py \
  --sv-file           data/my_svs.tsv \
  --gencode-file      annotation/gencode.v47.annotation.gtf.gz \
  --clingen-file      annotation/ClinGen_gene_curation_list.csv \
  --cosmic-file       annotation/cancer_gene_census.csv \
  --out-file          results/tables/fusion_final.txt \
  --html-index        results/report/index.html
```

Then open `results/report/index.html` in any browser.

---

## Input Formats

### Headerless legacy TSV

Seven tab-separated columns, no header row:

```
chrA  startA  endA  chrB  startB  endB  svtype
chr9  133713000  133713001  chr22  23632000  23632001  BND
```

`svtype` values: `DEL`, `DUP`, `INV`, `BND`

### MAVIS-style headered TSV

Must contain at minimum:

| Column | Description |
|---|---|
| `break1_chromosome` | Chromosome of breakpoint 1 |
| `break1_position_start` | Start of breakpoint 1 uncertainty interval |
| `break1_position_end` | End of breakpoint 1 uncertainty interval |
| `break2_chromosome` | Chromosome of breakpoint 2 |
| `break2_position_start` | Start of breakpoint 2 uncertainty interval |
| `break2_position_end` | End of breakpoint 2 uncertainty interval |
| `event_type` or `svtype` | SV type |
| `break1_orientation` | Breakpoint 1 strand orientation (`+` / `-`) — optional but recommended |
| `break2_orientation` | Breakpoint 2 strand orientation (`+` / `-`) — optional but recommended |

Uncertain breakpoints are resolved orientation-aware: `+` strand uses the lower coordinate, `-` strand uses the higher coordinate.

---

## Annotation Files

| File | Source | Purpose |
|---|---|---|
| GENCODE GTF (`.gtf.gz`) | [GENCODE](https://www.gencodegenes.org/human/) | Gene boundaries + CDS exons for reading-frame prediction. MANE Select and Ensembl canonical transcripts are preferred. |
| ClinGen dosage sensitivity | [ClinGen FTP](https://ftp.clinicalgenome.org/nightly/ClinGen_gene_curation_list_GRCh38.tsv) or browser CSV export | Haploinsufficiency (HI) and triplosensitivity (TS) scores used in ACMG scoring |
| COSMIC Cancer Gene Census | [COSMIC](https://cancer.sanger.ac.uk/census) (free registration) | Oncogene / tumour suppressor labels for `PM_ONCOGENE` and `BS_TSG_OOF` criteria |
| BED gene file | Custom or Ensembl/UCSC download | Alternative to GENCODE GTF; 4-column BED (chr, start, end, gene) with optional strand column 6 |

---

## Output Files

```
results/
 ├── tables/
 │   ├── fusion_final.txt              Main fusion summary (TSV, sorted by ACMG score)
 │   ├── variants_processed.tsv        All input SVs after loading and normalisation
 │   ├── breakpoint_annotations.tsv    Per-breakpoint gene overlap records
 │   ├── intragenic_events.tsv         Same-gene SVs excluded from main summary
 │   └── short_sv_events.tsv           Intrachromosomal SVs below --min-fusion-distance
 ├── report/
 │   └── index.html                    Self-contained interactive HTML report
 ├── svg/
 │   └── GENEA_GENEB_fusion.svg        One SVG diagram per fusion
 └── annotations/
     ├── fusion_breakpoints.bed        BED file for IGV / UCSC browser
     └── chimera_igv_session.xml       IGV session — File → Open Session
```

### fusion_final.txt columns

| Column | Description |
|---|---|
| `variant_ids` | Semicolon-separated variant IDs supporting this fusion |
| `geneA` / `geneB` | Gene symbols at each breakpoint |
| `break1_orientation` / `break2_orientation` | Breakpoint strand orientations |
| `strandA` / `strandB` | Gene strands from annotation |
| `strand_relation` | `SAME_STRAND` or `OPPOSITE_STRAND` |
| `svtypes` | SV type(s) supporting this fusion |
| `mechanism` | `DELETION_FUSION`, `DUPLICATION_FUSION`, `INVERSION_FUSION`, `TRANSLOCATION_FUSION`, `COMPLEX_SV_FUSION` |
| `fusion_class` | `CANONICAL`, `HIGH_CONF_CANONICAL`, `NON_CANONICAL`, `HIGH_CONF_NON_CANONICAL` |
| `reading_frame` | `IN_FRAME`, `OUT_OF_FRAME`, `UNKNOWN` |
| `support_score` | `0.5 × SV_count + 0.5 × unique_svtype_count` |
| `acmg_score` | Numeric ACMG score (see scoring table below) |
| `pathogenicity` | `Pathogenic`, `Likely_Pathogenic`, `Uncertain_Significance`, `Likely_Benign`, `Benign` |
| `acmg_codes` | Semicolon-separated ACMG criteria codes applied |
| `pathogenicity_evidence` | Human-readable summary of evidence |
| `reciprocal_pair_id` | Links reciprocal fusion pairs sharing a variant ID |

---

## ACMG Scoring

Chimera applies an ACMG/AMP-inspired point-based scoring system. Criteria are applied additively; the total score determines the pathogenicity call.

### Pathogenic criteria

| Code | Score | Condition |
|---|---|---|
| `PVS1` | +8 | Both genes haploinsufficient (ClinGen HI ≥ 2), canonical junction, in-frame |
| `PM_INFRAME` | +2 | In-frame junction confirmed (additive with HI criteria) |
| `PS1` | +6 | Both genes haploinsufficient, canonical, frame unknown or out-of-frame |
| `PM1` | +2 | One gene haploinsufficient, canonical |
| `PM_INFRAME_NOHI` | +4 | Canonical + in-frame, no ClinGen dosage data |
| `PM_ONCOGENE` | +2 | Known oncogene (COSMIC CGC), canonical fusion |
| `PP3` | +1 | Canonical junction, reading frame undetermined |
| `PP_WEAK_HI` | +1 | HI score = 1 (limited evidence) |

### Benign criteria

| Code | Score | Condition |
|---|---|---|
| `BP_OOF` | -2 | Confirmed out-of-frame junction |
| `BP4` | -1 | Both genes dosage-insensitive |
| `BP7` | -1 | Same gene at both breakpoints (intragenic artifact) |
| `BP_AR` | -1 | Gene is autosomal recessive only (ClinGen HI score 30) |
| `BS1` | -4 | Duplication; neither gene triplosensitive |
| `BS_TSG_OOF` | -2 | Tumour suppressor gene + out-of-frame junction |

Non-canonical junctions receive an additional **−3 penalty** regardless of other criteria.

### Score → classification

| Score | Classification |
|---|---|
| ≥ 10 | Pathogenic |
| 6 – 9 | Likely Pathogenic |
| −5 to 5 | Uncertain Significance |
| −6 to −9 | Likely Benign |
| ≤ −10 | Benign |

---

## All Command-Line Options

```
python chimera.py [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `--sv-file` | `Book2.clean.tsv` | Input SV file (headerless TSV or MAVIS-style headered TSV) |
| `--gene-file` | `annotation/genes.bed` | BED gene annotation (used if `--gencode-file` is absent) |
| `--gencode-file` | `annotation/gencode.v47.annotation.gtf.gz` | GENCODE GTF/GFF3 (`.gz` supported). Takes precedence over `--gene-file` |
| `--out-file` | `results/tables/fusion_final.txt` | Main fusion summary output |
| `--variants-out` | `results/tables/variants_processed.tsv` | Processed SV table |
| `--breakpoints-out` | `results/tables/breakpoint_annotations.tsv` | Per-breakpoint gene annotations |
| `--svg-dir` | `results/svg/` | Directory for SVG fusion diagrams |
| `--max-svgs` | `500` | Maximum number of SVG diagrams to generate |
| `--svg-width` | `1200` | SVG width in pixels |
| `--clingen-file` | `annotation/tableExport-2.csv` | ClinGen dosage sensitivity CSV/TSV |
| `--cosmic-file` | `annotation/cancer_gene_census.csv` | COSMIC Cancer Gene Census CSV/TSV |
| `--min-sv-score` | `0.0` | Minimum `support_score` to include a fusion in output |
| `--min-fusion-distance` | `10000` | Minimum bp distance for intrachromosomal fusions; closer events go to `--short-sv-out` |
| `--intragenic-out` | `results/tables/intragenic_events.tsv` | Output for same-gene breakpoint events |
| `--short-sv-out` | `results/tables/short_sv_events.tsv` | Output for near-distance intrachromosomal events |
| `--html-index` | `results/report/index.html` | Self-contained HTML report |
| `--cohort-report` | `results/report/cohort_summary.html` | Standalone cohort summary (legacy; merged into `--html-index`) |
| `--bed-out` | `results/annotations/fusion_breakpoints.bed` | BED file for genome browser loading |
| `--igv-session` | `results/annotations/chimera_igv_session.xml` | IGV session XML |
| `--genome` | `hg38` | Reference genome label for IGV session |
| `--debug-genes` | _(none)_ | One or more gene symbols for CDS lookup diagnostics, e.g. `--debug-genes ABL1 BCR` |

---

## Project Layout

```
chimera/
 ├── chimera.py                  Pipeline script
 ├── README.md
 ├── annotation/
 │   ├── gencode.v47.annotation.gtf.gz
 │   ├── ClinGen_gene_curation_list.csv
 │   └── cancer_gene_census.csv
 ├── data/
 │   └── my_svs.tsv
 └── results/                    Created automatically on first run
```

---

## How It Works

```
SV file
   │
   ▼
Load & normalise SVs ──────────────────────────────────────────────┐
   │  (headerless or MAVIS format; uncertain breakpoints resolved)  │
   ▼                                                                │
Overlap breakpoints with gene index                                 │
   │  (exact interval overlap; no window expansion)                 │
   ▼                                                                │
Build fusion edge graph                                             │
   │  • same-gene events → intragenic_events.tsv                   │
   │  • too-close pairs  → short_sv_events.tsv                     │
   ▼                                                                │
Classify each fusion                                                │
   │  • canonical/non-canonical via strand × SV-type matrix        │
   │  • reading frame via GENCODE CDS in superTranscript space      │
   │  • ACMG score via ClinGen HI/TS + COSMIC oncogene/TSG         │
   ▼                                                                │
Write outputs                                                       │
   │  • fusion_final.txt (TSV)                                      │
   │  • SVG diagrams                                                │
   │  • index.html (self-contained report)                          │
   │  • fusion_breakpoints.bed + chimera_igv_session.xml            │
   └───────────────────────────────────────────────────────────────►
```

---

## HTML Report Tabs

| Tab | Contents |
|---|---|
| **Cohort Summary** | KPI cards, pathogenicity pie charts, ACMG score histogram, SV mechanism bars, top genes table, priority fusions table |
| **Fusion Annotations** | One collapsible card per fusion with inline SVG diagram, metadata table, ACMG breakdown bar, compare mode |
| **Visualisations** | Circos-style genome breakpoint map, ACMG rank plot, force-directed fusion network |
| **Raw Data** | Searchable/sortable/downloadable tables for all five output TSVs |
| **Files & Tools** | Output file listing with sizes, IGV and UCSC browser integration instructions, ACMG criteria glossary |

---

## Genome Browser Integration

**IGV Desktop**
```
File → Open Session → chimera_igv_session.xml
```
IGV opens at the top-priority fusion locus with breakpoints pre-loaded and colour-coded by pathogenicity.

**UCSC Browser**
Host `fusion_breakpoints.bed` on a web-accessible server, then go to *My Data → Custom Tracks* and paste the URL.

---

## Citation

If you use Chimera in your research, please cite this repository and acknowledge the following data sources:

- Rehm et al. (2015) ClinGen — The Clinical Genome Resource. *NEJM* 372:2235–2242
- Sondka et al. (2018) The COSMIC Cancer Gene Census. *Nature Reviews Cancer* 18:696–705
- GENCODE: Frankish et al. (2021) *Nucleic Acids Research* 49:D916–D923

---

## License

MIT
