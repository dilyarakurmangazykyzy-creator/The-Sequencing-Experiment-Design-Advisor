# Scientific sources and interpretation rules

Project 15: **The Sequencing Experiment Design Advisor**. Sources checked on 2026-10-06. This document records evidence for the design and the limits of the implemented bacterial pilot. It distinguishes measured results, mathematical predictions, published results and editable planning assumptions.

## The biological question

Given an organism, a biological objective and a budget, which sequencing design is feasible, and what evidence supports its read length and depth? Genome size determines the amount of sequence needed for a nominal DNA depth; ploidy, allele fraction, repeats, sample complexity, library biases and the endpoint determine whether that depth is useful.

The empirical component of this project is a bacterial DNA pilot. Human variants, RNA expression and finished assemblies require separate endpoint-specific validation. A bacterial coverage curve cannot establish a universal depth threshold for all applications.

| Objective | Why depth matters | Why depth alone is insufficient | Suitable saturation endpoint |
|---|---|---|---|
| Small variants | More independent observations increase support for an allele | Reference divergence, repeats, systematic errors, contamination and caller thresholds can dominate | Precision, recall and F1 against an independent truth set; stratify difficult regions |
| Gene expression | Deeper libraries sample low-abundance transcripts more often | RNA is sampled by abundance, not uniformly across genome bases; biological variability requires biological replicates | Detected genes, stable count estimates and differential-expression power across depths and replicates |
| Assembly | More overlaps can increase contiguity and consensus support | Repeats longer than the informative molecule span may remain unresolved at arbitrarily high short-read depth | Genome fraction, misassemblies, NGA50 and sequence errors; N50 alone is insufficient |

Primary empirical evidence: Liu et al. analysed depth versus replication in RNA-seq and found diminishing returns from additional depth in their MCF7 experiment; their threshold is specific to those data. [RNA-seq differential expression studies: more sequence or more replication?](https://pmc.ncbi.nlm.nih.gov/articles/PMC3904521/). Schurch et al. directly evaluated biological replicate counts in a large yeast experiment. [How many biological replicates are needed in an RNA-seq experiment?](https://pmc.ncbi.nlm.nih.gov/articles/PMC4878611/).

Assembly evidence: Wick et al. compared bacterial assemblies using short and long reads and showed that resolving repeats changes the outcome even when short-read coverage is abundant. They also documented genuine subculture differences between MG1655 reads and its historical reference. [Unicycler](https://doi.org/10.1371/journal.pcbi.1005595). Bresler et al. derived reconstruction limits depending jointly on read length, coverage and repeat structure. [Optimal assembly for high throughput shotgun sequencing](https://doi.org/10.1186/1471-2105-14-S5-S18).

Long-read variant evidence: Harvey et al. performed downsampling across long-read technologies, depths and read/assembly-based methods. This supports testing the actual endpoint and analysis method rather than transferring a generic coverage threshold. Use the published article rather than treating its earlier PMC preprint as the final paper. [Genome Research published article](https://doi.org/10.1101/gr.278070.123).

## Definitions used in the analysis

Let G be the total length of the selected reference, B the number of sequenced bases, and d_i the accepted aligned depth at reference position i.

| Metric | Definition | Interpretation and caveat |
|---|---|---|
| Nominal/raw depth | B/G, in × | Includes unmapped, low-quality and duplicate sequence; not usable local depth |
| Pilot accepted-molecule nominal depth | Original raw bases of molecules accepted by QC / G in the SAM route | Distinct from all acquired raw depth. SAM intervals describe the original aligned input, so audited trimming does not shorten those intervals |
| Mean aligned depth | sum(d_i)/G, including zero-depth positions | Requires a specified CIGAR, flag, base-quality and mapping-quality policy |
| Breadth at k | count(d_i >= k)/G | Fraction of the reference with at least k accepted observations |
| Mapping fraction | Mapped primary read records / all primary read records | State whether the unit is a read or a paired molecule; exclude repeated counting of secondary/supplementary records |
| Alignment identity | Matches / alignment columns under the stated match/mismatch/indel convention | Requires a real alignment; matching anchors alone do not provide nucleotide identity |
| Mapping quality | -10 log10(probability the placement is incorrect), when available | Aligner estimate, not measured accuracy; SAM value 255 means unavailable |
| Q20/Q30 base fraction | Fraction of bases with Phred >=20 / >=30 | Quality labels are model predictions; empirical error requires comparison to trusted sequence |
| Expected error fraction | mean over bases of 10^(-Q/10) | Average probabilities, not 10^(-mean Q/10) |
| Read N50 | Length L such that reads of length >=L contain at least half the bases | A length statistic, not accuracy or an assembly quality statistic |

The official SAM specification defines Phred scale, mapping quality, alignment flags, CIGAR and coordinate conventions. [SAM/BAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf). Samtools explicitly distinguishes covered-base percentage from mean depth. [samtools coverage](https://www.htslib.org/doc/samtools-coverage.html).

For standard coverage, zero-depth positions must remain in the denominator. Specify paired-end overlap treatment, deletion treatment, duplicate flags and thresholds. Note that `samtools depth -q` filters **base** quality and `-Q` filters **mapping** quality; those letters are reversed for `samtools coverage`. [samtools depth](https://www.htslib.org/doc/samtools-depth.html).

### The mapping-rate premise in assignment task 6

The assignment asks how mapping rate and alignment quality degrade as coverage falls. That is a hypothesis to test. When each read is aligned independently against a fixed reference with fixed parameters, a uniformly retained subset should have approximately the same expected mapping fraction. Its sampling uncertainty grows as the number of reads falls. Coverage breadth, allele support and the reliability of an answer can decline strongly while mapping fraction remains high.

For independent sampling, if p is the full-pool mapping fraction and n reads are retained, the standard error of the observed fraction is approximately sqrt(p(1-p)/n). Paired molecules violate independence between mates, so molecule-level resampling is preferable for uncertainty estimates. This is a mathematical explanation, not a result attributed to an external paper.

Different mapping fractions after downsampling can also come from biased sampling, changed filters, different library composition, or algorithms that share information across reads. A nearly flat curve is an informative negative result. Do not alter a curve to match the wording of the brief.

### What the dependency-free mapper measures

The fallback in `pipeline.py` is an educational seed/anchor mapper. Short reads use candidate mismatch/banded-edit verification with approximate placement intervals; long reads use conservative exact-anchor-supported intervals. Its mapping fraction is acceptance by that algorithm, and its breadth/depth are **proxy coverage**. It cannot supply calibrated MAPQ, general CIGAR identity, variant calls, SV sensitivity or assembly correctness. Exact anchors preferentially miss error-rich reads, making platform comparisons partly algorithm-dependent.

The completed standard pilot uses genuine minimap2 2.22 WebAssembly SIMD, `-ax sr` for paired short DNA reads and `-ax map-ont` for ONT reads, together with primary SAM output and CIGAR-aware processing. `map-ont` is an ONT error-profile preset available in this version; the archived 2025 experiment's chemistry is not established from its instrument name alone. Recent high-accuracy ONT has a different preset (`lr:hq`); old PacBio CLR differs from HiFi. Record the tool version and all presets. [minimap2 official README](https://github.com/lh3/minimap2), [minimap2 manual](https://github.com/lh3/minimap2/blob/master/minimap2.1).

## Downsampling and saturation design

Downsampling must retain both Illumina mates using one molecule identifier and one deterministic random decision. Within a seed, nested subsets reduce arbitrary comparisons; multiple seeds estimate sampling variability. Report raw/filtered reads and bases, realized rather than merely requested depth, retained fraction and uncertainty at every point.

An ordered prefix obtained by streaming a large file is an **acquisition cap**, not a random sample of the full run. It may contain a non-representative sequencing-time, tile or quality distribution. Randomly downsample within that acquired pool and explicitly limit conclusions to it. A full-file streaming reservoir/hash sample is a better extension when transfer budget allows.

Define saturation before choosing a result. A transparent pilot rule is: choose the lowest observed depth whose mean measured breadth passes the required target and whose subsequent increments change that endpoint by less than a fixed epsilon. State target, epsilon, seed count and how many later points are required. If no tested depth meets the rule, return **not reached**. The highest tested value alone is not proof of a plateau.

For an accuracy claim, use independent truth rather than the full-depth calls as unquestioned truth. Full-depth calls can instead be labelled a concordance baseline. Hold out a run/sample for testing when fitting a saturation or recommendation model; repeated subsamples from one pool are technical resampling, not biological replication.

## Predicted detection support

The analytical advisor can estimate the probability of enough supporting observations. If accepted local depth is fixed at n and the true allele fraction is f, K~Binomial(n,f). When depth is represented as a Poisson mean D_eff, K~Poisson(D_eff*f). For a requirement of k alternative observations:

`P(K >= k) = 1 - exp(-D_eff*f) * sum((D_eff*f)^j / j!, j=0..k-1)`.

For k=3 and a 95% support target, the required mean supporting count is approximately 6.296. Thus the idealized fraction threshold is approximately `f95 = 6.296/D_eff`; if it exceeds 1, the target is infeasible. For k=1 it is approximately `2.996/D_eff`. These values follow directly from the probability model.

This is a **sampling support probability**, not a truth-validated limit of detection. It excludes systematic errors, mapping bias, PCR duplicate dependence, strand balance, local depth variation and false-positive control. A caller can miss an allele despite sufficient support, or call an artifact supported by many reads. Clinical or biological detection accuracy needs its own truth-set benchmark.

A clonal haploid bacterial SNV generally has f near 1; f=0.5 represents a mixed population/allele mixture, not a diploid bacterial heterozygote. For a diploid human germline heterozygote f near 0.5 is an ideal starting assumption. Metagenome allele support must additionally account for the organism's relative abundance.

## Archives, identifiers and dataset evidence

Use explicit joins: run -> experiment -> sample/BioSample -> study/BioProject -> publication/reference. A strain name alone does not prove that two technologies sequenced the same aliquot. A shared BioSample is stronger evidence, although laboratory handling can still differ. For cross-archive disagreement, preserve both raw records and report which field is used and why.

| Selected dataset | Evidence | Limits and acquisition considerations |
|---|---|---|
| MG1655 Illumina `SRR1030394` | BioSample `SAMN02401358`, BioProject `PRJNA227741`, Illumina MiSeq, paired WGS. [Published benchmark identifying this run](https://academic.oup.com/view-large/409613964) | ENA first-public 2014-11-27. Complete source contains 2,720,956 spots/pairs and 1,322,136,908 bases. Acquired prefix is smaller and must retain both mates |
| MG1655 ONT `SRR33602302` | Experiment `SRX28830916`, BioSample `SAMN48541625`, BioProject `PRJNA1222438`, SRA study `SRP563276`, GridION, single reads. [NCBI primary experiment record](https://www.ncbi.nlm.nih.gov/sra/SRX28830916%5Baccn%5D) | Complete source reports 115,315 reads and 410,399,976 bases. Library Strategy says **WGA**, while Design says **WGS**. Actual amplification status is therefore uncertain from these metadata alone |
| Reference `U00096.3` | GenBank counterpart of MG1655 RefSeq `NC_000913.3`, retrieved through ENA browser API | A strain reference does not establish exact biological truth for different subcultures. Record FASTA checksum and preserve accession version |

The actual acquired count, retrieval date, checksum and source URL are recorded in `data/raw/manifest.json` and the saved provenance copies. Complete-source counts above are **archive metadata**, not measured acquired counts. The short/long runs share a named strain but different BioSamples, projects and years. They cannot isolate the platform effect.

Concrete metadata issues: the ONT NCBI record labels Design WGS and Strategy WGA, so preserve both fields and describe the conflict. NCBI reports publication 2025-05-23 while ENA reports first-public 2025-05-24, so report the archive field actually used instead of silently reconciling the one-day difference. ENA Illumina `base_count/read_count` is approximately 486 bp in this run and represents the combined spot/pair span, while individual mates average roughly 243 bp. Check actual FASTQ structure before interpreting a generic count as individual reads.

### Scalable metadata retrieval

The ENA Portal API supports `result=read_run`, field selection, query filters and TSV output. Its default search ceiling is 100,000 rows; `limit=0` requests all matching records, including sets above one million. ENA advises no more than ten simultaneous searches. A reproducible query plus streaming parsing provides a scalable retrieval method; only the actual acquired row count may be described as completed. [ENA advanced programmatic search](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/advanced-search.html).

NCBI provides the SRA metadata in BigQuery and Athena. Its schema includes accession, platform, instrument, assay type, library layout, BioSample, BioProject and release date. [SRA cloud metadata schema](https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-based-metadata-table/), [cloud querying overview](https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-based-examples/).

The official BigQuery dataset is `nih-sra-datastore.sra.metadata`. Full-archive grouping can aggregate without downloading all raw rows, but needs the user's cloud project/account. Save SQL, job metadata, bytes processed, returned groups, date and table schema. [NCBI BigQuery setup](https://www.ncbi.nlm.nih.gov/sra/docs/sra-bigquery/).

Metadata year must be labelled **first-public/release year**, unless actual sequencing date is available. Library strategy is an application proxy, not proof of study intent or successful outcomes. Missing depth cannot be replaced by zero. `bases / spots` can be mean **spot** length; for paired data it must not automatically be called read length. Inspect read structure and report assumptions before converting. [SRA Handbook](https://www.ncbi.nlm.nih.gov/books/NBK47528/pdf/Bookshelf_NBK47528.pdf).

Publication-linked selection has survivorship and reporting bias. Archives describe what was submitted, not necessarily what succeeded, nor the price paid. Include missingness and inconsistent-field counts, deduplicate by run and identify archived mirrors rather than double-counting them as independent observations.

## Instrument and cost evidence

Current advertised throughput is conditional on chemistry, flow-cell configuration and sample/library quality. It is a capability envelope, not a guaranteed usable yield. Current prices cannot be attributed to archived experiments from a different year.

| Source | Evidence it supports | Scope |
|---|---|---|
| [Illumina NovaSeq X specifications](https://www.illumina.com/systems/sequencing-platforms/novaseq-x-plus/specifications.html) | Single 25B flow cell at PE150: approximately 8–10.5 Tb; dual-cell X Plus doubles single-cell output | Vendor specification; upper range not guaranteed |
| [PacBio HiFi grant support and specifications](https://www.pacb.com/hifi-grant-support/) | Revio + SPRQ-Nx advertises 120–480 Gb per run depending on acquisitions; HiFi read accuracy and current platform configurations | Vendor claims; keep CLR/HiFi distinct and name chemistry |
| [Oxford Nanopore chemistry technical document](https://nanoporetech.com/document/chemistry-technical-document) | Read length and accuracy depend on chemistry, extraction, size selection and basecalling | Use exact experimental conditions; instrument and year alone do not identify the archived pilot's chemistry |
| [Oxford Nanopore US store](https://store.nanoporetech.com/us/flow-cells.html) | Single FLO-MIN114 DNA flow cell: US$840 as displayed on the check date | Flow-cell purchase only; prep, device, tax and labor excluded |
| [Northwestern NUSeq FY26 service prices](https://www.cgm.northwestern.edu/cores/nuseq/pricing.html) | Independent public service benchmarks in USD: external NovaSeq X 10B PE150 lane $1,873.75 for 1.0–1.2B read pairs; 25B PE150 lane $3,750 for 2.5–3.0B pairs; MinION $1,250/cell up to50Gb; ONT prep $437.50/sample; WGS prep <24 samples $240/sample | External academic service rates; capacity and minimum billable unit matter; quotation needed for actual projects |
| [NHGRI sequencing cost accounting](https://www.genome.gov/about-genomics/fact-sheets/DNA-Sequencing-Costs-Data) | Separates production cost metrics and explains included/excluded cost scope | Historical benchmark; not a platform-specific local quote |

A cost model should add preparation/QC + sequencing consumables + allocated instrument time + analysis compute + storage/retention + specified contingency. A service price may already include consumables and instrument time; do not add those a second time. Shared-pool and dedicated-run estimates must be separate. Billable runs are discrete: `ceil(required_raw_Gb/usable_Gb_per_run)` rather than a fraction of a cell for an unshared project.

**Independent validation example:** 1.0–1.2 billion PE150 pairs contain 300–360 Gb. The NUSeq external 10B lane price therefore implies about $5.205–$6.246 per raw Gb for the sequencing service alone. The 25B PE150 lane contains750–900Gb and implies $4.167–$5.000/Gb; it matches the model's selected 25B configuration more closely. Compare the model's sequencing-plus-instrument prediction with the corresponding interval and report the relative difference, while preserving the facility's whole-lane minimum. This is a derived service benchmark, not validation of preparation, usable yield, compute or total cost. Matching a vendor number used as the model's own input is an arithmetic consistency check, not independent validation.

Reagent list prices, independent service prices and editable assumptions must retain separate provenance. Record currency, region, date, instrument/configuration, minimum order, included items and yield assumption. Sensitivity analysis should vary usable yield, sharing, DNA/RNA library costs and compute. A fixed ±30% planning band is a scenario range, not a statistically fitted confidence interval.

## Claims permitted by the current evidence

- The executed pilot uses 400,000 real archived Illumina pairs and 60,000 real archived ONT reads; the small included fixture verifies the workflow and does not reproduce the deep-pool plateau.
- Genuine minimap2 primary SAM/CIGAR intervals support measured reference mapping and breadth with the stated MAPQ/flag/paired-overlap policy. The separate educational fallback retains its proxy label.
- The acquired pools yield operational reference-breadth plateaus at 14.121× Illumina and 15.810× ONT under the declared rule. These thresholds are not variant recall or biological detection limits.
- `results/saturation.json` measures mapping fraction, available MAPQ and NM/aligned-query edit-distance proxy at nine depths and three technical seeds. Stable mapping fraction is an informative negative finding. At full depth, NM proxy is approximately 0.9523% for Illumina and 0.6076% for ONT; different samples/libraries/years prevent a causal platform accuracy ranking.
- Streaming ENA metadata aggregation processed 999,988 actual run records across 12 platforms. The unsorted cap selects only 2022–2025 records from a 2015–2025 query, so it cannot establish full-archive prevalence or an unbiased temporal trend. ENA and SRA are archive mirrors rather than independent replicated observations.
- The recommender returns a justified planning design, explicit cost components and analytical observation-support probability. Independent NUSeq comparisons cover selected sequencing/preparation components and reveal scope/price mismatches; preparation, compute and project totals remain conditional.
- RNA power, human variant accuracy and complete repeat assemblies need their corresponding endpoint-specific analyses before being described as validated. This is a limitation on those claims, not a statement that every such additional experiment is mandated for the bacterial breadth pilot.

The exact values, scope and provenance are in `results/summary.json`, `results/saturation.json`, `results/metadata_summary.json`, `results/acquisition_manifest.json` and `results/metadata_provenance.json`. The executed accepted-molecule/raw-alignment route is explained in `docs/ALIGNMENT_METHOD.md`; the prepared native fastp/minimap2/SAMtools workflow has not been executed here.

Do not promise a grade. The course handbook's rubric totals 60 points; excellent files do not replace a genuine two-person contribution record or understanding at the defence.
