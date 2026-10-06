# The Sequencing Experiment Design Advisor

AITU Introduction to Bioinformatics, Project 15. Team: **Zhaulybaeva Gaukhar and Kurmangazykyzy Dilyara**.

The advisor selects a platform, read length and planned depth under a budget. It shows the equations, editable costs, alternatives and limits of its evidence. Real bacterial Illumina and Nanopore reads support reference-coverage saturation experiments. Statistical detection support and other application recommendations are explicitly labeled planning models.

Submission documents: [report with team responsibilities and Appendix A: AI assistance](docs/report.pdf), [12-slide presentation with bilingual speaker notes](docs/presentation.pptx), and [contribution record](docs/CONTRIBUTIONS.md).

Each member can run her own test launcher (`CHECK_GAUKHAR.cmd` / `CHECK_DILYARA.cmd`) from the complete project and document her actual results. The launchers save execution records in `results/reviews/`. Repository uploads and independent technical checks are recorded separately; AI assistance is disclosed in Appendix A of the report.

## Launch in under ten steps

1. Install Python 3.12 or newer.
2. Download/clone this repository and open a terminal in its folder.
3. Run `python app.py` (Windows: double-click `START_ADVISOR.cmd`).
4. Open http://127.0.0.1:8000.
5. Pick an example or set organism, biological question, sample count and budget.
6. Inspect the feasible design, assumptions, measured bacterial evidence and alternatives.
7. Export the design JSON for your report.

The advisor has no third-party dependencies and no cloud API keys. No user-supplied biological data are transmitted. A Node.js installation is needed only for the WebAssembly alignment experiment.

## Delivered measured results

The completed bounded cohort contains **999,988 ENA run records** across 12 platforms. Streaming aggregation took 50.40 seconds with 36.82 MiB peak process memory. This is an unsorted capped query cohort, not the entire archive.

The real pilot acquired **400,000 Illumina pairs and 60,000 Nanopore reads**. QC retained 384,691 pairs and 50,175 long reads. Primary mapping rates in those accepted pools were 99.6020% and 99.9821%. MAPQ>=20 reference breadth at >=1x was 99.1379% and 99.9378%. The declared reference-breadth plateau was 14.12x and 15.81x accepted-molecule raw nominal depth. These values do not establish a universal sequencing requirement or a platform accuracy ranking.

All 30 unit tests passed, and the documented one-command small real-read verification completed using checksum-pinned minimap2 SIMD assets. `results/` includes compact QC, saturation, metadata, cost examples and resource/provenance summaries.

## Reproduce the real-data experiment

Install Node.js 18 or newer. The included small **real archive read subset** verifies retrieval-independent alignment and analysis in one command. The first run downloads checksum-pinned minimap2 2.22 WebAssembly assets (~1 MB).

```text
python run_experiment.py
```

This writes `results/quick/`. Its tiny read set checks that the pipeline works; it does not reproduce the large pilot's plateau. To reproduce the headline pilot from raw archived input:

```text
python -m pip install -r requirements.txt
python run_experiment.py --pilot
```

The second command downloads 400,000 Illumina pairs and 60,000 Nanopore reads (~350 MB compressed retained input), a 4.64-Mb reference and metadata, aligns with `sr`/`map-ont`, performs three seeded nested downsampling replicates, and writes `results/summary.json` and `results/saturation.json`. Prefix acquisition is disclosed in the manifest. Subsampling within those prefixes is deterministic; it is not a random sample of the complete source runs. Fresh archive queries can change, so preserve the query, retrieval timestamp, SHA-256 and measured row counts.

For metadata at scale (streaming, bounded cohort; `limit=0` requests all matching rows):

```text
python scripts/metadata_mine.py --fetch --input data/raw/metadata_million.tsv --out results --query "first_public>=2015-01-01 AND first_public<=2025-12-31" --limit 1000000
```

The delivered `metadata_summary.json` records the actual completed query, not a promised archive size. This cohort is unsorted and cannot establish representative platform market shares or experimental success rates.

The Dockerfile specifies a reproducible Python environment for CLI analysis: `docker build -t sequencing-advisor .` then `docker run --rm sequencing-advisor python run_experiment.py`. Docker was unavailable on the development host, so the image was not executed. The native Linux Snakemake route realigns each fraction with fastp/minimap2/samtools and records resources:

```text
snakemake --snakefile workflow/Snakefile --cores 4 --use-conda
```

Install Snakemake in a Linux/WSL or institutional environment before that optional command. The executed route is Python+Node minimap2 WebAssembly; the native workflow is supplied and not claimed as an executed result.

## Data and interpretation

| Component | Accession | What was acquired |
|---|---|---|
| Illumina MiSeq WGS | SRR1030394 / SAMN02401358 / PRJNA227741 | First 400,000 archived read pairs |
| Nanopore GridION | SRR33602302 / SAMN48541625 / PRJNA1222438 | First 60,000 archived reads |
| MG1655 reference | U00096.3 (GenBank), NC_000913.3 counterpart | Complete reference sequence |

Both runs name MG1655, but they are different samples/studies/years. ENA labels the Nanopore library WGA while NCBI's Design text says WGS: amplification status is uncertain. The comparison does not isolate a platform effect. The paper linked from the Nanopore BioProject studies retron biology, not a platform benchmark.

Mapping rate counts primary mapped reads. Coverage uses primary aligned M/=X bases at MAPQ >=20, excluding flagged duplicate/QC-fail reads and MAPQ255. Nominal depth uses acquired raw bases divided by reference length. Paired overlapping bases count twice; this is read depth, not unique-molecule depth. QC-selected molecules retain their original raw-read alignments in the executed route. The native workflow realigns trimmed reads.

Saturation means >=95% reference breadth at >=1x and two adjacent gains below one percentage point. This endpoint does not validate variant recall, expression accuracy or assembly completeness. Decreasing depth removes reference support; per-read mapping rate need not decrease. The Poisson detection calculation predicts alternate-read support, not a truth-validated clinical limit of detection.

## Checks and files

```text
python -m unittest discover -s tests -v
python advisor.py --organism bacteria --application variant --budget 1200 --samples 6
```

`pipeline.py` handles FASTQ QC, paired identities, SAM/CIGAR coverage and nested experiments. `advisor.py` contains the transparent recommendation model. `config/costs.json` contains source-backed components and clearly labeled assumptions. `scripts/metadata_mine.py` streams metadata. `docs/` contains the report, presentation, contribution record, scientific sources and alignment methods.

Large reads and alignments stay outside Git. Include the small real `data/example/` fixture, compact result JSON/CSV, source code, report and slides. Both partners must make actual technical contributions and explain all components; Appendix A of the report describes the AI assistance used.
compact result JSON/CSV, source code, report and slides. Both partners must make actual technical contributions and explain all components; AI assistance is disclosed in `docs/AI_DISCLOSURE.md`.
