# Completion audit against the course handbook

This audit is a submission checklist, not a predicted grade. The supplied handbook assigns **60 points** to the course project and requires real archived data, an 8–12 page report, a reproducibility bundle, a two-person contribution statement and AI disclosure. The examiner decides the score. There is no defensible guarantee of “100/100”.

Status meanings: **Implemented** = a capability is present in the source; **Measured** = an output from an actual run is present and checked; **Partial** = useful evidence exists but the stated full task is broader; **Required from team** = names, authentic contributions, submission or defence cannot be completed by an assistant. An implementation or a prepared command is not a completed experiment.

This review checks the actual final pilot outputs. It distinguishes that executed route from the separately prepared native fastp/minimap2/SAMtools workflow. The final report and presentation must preserve the same empirical scope as those outputs.

## Checked empirical evidence

- `results/summary.json`: 400,000 Illumina paired molecules and 60,000 ONT reads against the 4,641,652-bp MG1655 reference. QC retains 384,691 pairs and 50,175 ONT reads. Minimap2 2.22 WebAssembly SIMD, `sr` and `map-ont`, produces genuine primary SAM/CIGAR alignments.
- Accepted-molecule nominal depth is 40.346× and 45.047×; mapped primary fractions are 99.6020% and 99.9821%; reference breadth at ≥1× with MAPQ≥20 is 99.1379% and 99.9378%, respectively. The nominal numerator uses original raw lengths of QC-accepted molecules. The alignments were not regenerated after audited trimming.
- `results/saturation.json`: nine nested fractions, three technical seeds. The declared breadth plateau is 14.121× for Illumina and 15.810× for ONT. Mapping fraction, available MAPQ and NM/aligned-query edit-distance proxy are measured at each depth. These are not variant-caller detection accuracy.
- `results/metadata_summary.json` and `results/metadata_designs.tsv`: 999,988 actual ENA run records processed from a capped, unsorted 2015–2025 query. Observed records span 2022–2025 and 12 platforms; three displayed DNA platforms account for 970,678 records. This cap does not establish full-archive frequencies.
- `results/acquisition_manifest.json` and `results/metadata_provenance.json`: source URLs, query/date and local checksums. Raw prefixes are genuine ordered acquisitions; full-source archive checksums were not verified. Alignment version/preset/resource logs are produced locally and summarized in the report.
- The main QC/coverage/downsampling pass took 244.74 s and observed 748.32 MiB peak process memory; this excludes acquisition and alignment. Outer JavaScript alignment elapsed times are 181.38 s for Illumina and 221.94 s for ONT; internal minimap2 logs report 181.35 s and 221.93 s. Completion RSS is not an alignment peak-memory measurement. The metadata aggregation took 50.40 s, processed about 19,841 records/s and observed 36.82 MiB peak process memory; `results/metadata_scale.json` records fresh-process scaling checks.

## Eight project-specific tasks

| Task from Project 15 | Current evidence and intended file | Status / limitation | What would strengthen completion |
|---|---|---|---|
| 1. Explain design limits and different depth requirements for variants, expression and assembly | `docs/SCIENTIFIC_SOURCES.md`; report biological framing; advisor explanations | Implemented scientific framing. The empirical pilot is bacterial DNA | Explain ploidy, allele fraction, library complexity, uneven expression, repeat span and independent biological replicates at defence |
| 2. Mine SRA metadata at scale by application and year | `scripts/metadata_mine.py`; `results/metadata_summary.json`, `metadata_designs.tsv`, `metadata_provenance.json` | Measured 999,988-run streaming aggregation; partial full-archive scope. Library strategy is an application proxy, and year is archive release year | Unlimited/partitioned ENA or cloud aggregation with exact query/schema, deduplication, missingness and study/publication curation. Do not infer population trends from the cap |
| 3. Model output/read length/error/cost envelopes | `config/costs.json`, `advisor.py`, source inventory | Implemented planning model; specifications and prices are conditional on configuration/date; archived chemistry differs from current instruments | Validate configurations, maintain units, separate raw from usable yield, and include minimum billable run and shared-pool behavior |
| 4. Downsample deep short reads and find saturation | `pipeline.py`; `results/summary.json`, `saturation.json`; report/deck breadth curves | Measured 40.346× accepted raw nominal pool, genuine CIGAR coverage, three seeds and declared 14.121× breadth plateau. Limited to the acquired prefix and breadth endpoint | A complete/randomized-source acquisition improves representativeness. Independent truth is necessary for additional variant/expression/assembly accuracy claims, rather than for this stated breadth endpoint |
| 5. Repeat for long reads and demonstrate where they are necessary | Same results; separate ONT filters; literature and advisor repeat-span rules | Measured 45.047× pool and 15.810× breadth plateau. Literature/model explain read-span necessity; empirical repeat assembly/SV proof remains partial | Matched samples and a repeat/assembly or SV endpoint would demonstrate the structural advantage directly. Current samples cannot isolate a platform effect |
| 6. Measure alignment quality/mapping rate versus depth | Genuine minimap2 SAM; `saturation.json` per-depth mapping/MAPQ/NM means and seed ranges; `docs/ALIGNMENT_METHOD.md` | Measured. Approximately stable mapping fraction is an informative negative result. MAPQ is uncalibrated and NM is an edit-distance proxy | Use independent truth for placement/base-error accuracy. Retain the actual flag/MAPQ/CIGAR and overlapping-pair policy; do not fabricate a mapping decline |
| 7. Combine consumables, instrument time, cloud cost and validate | `config/costs.json`, `advisor.py`; `results/recommendation_*.json` cost breakdown and `cost_validation` | Implemented and compared with independent NUSeq sequencing/preparation services. Vendor consistency and facility comparisons have different scopes; total cost is not independently validated | Obtain a local matching-scope quote and measured compute/storage service cost. Include preparation-price mismatch and shared-pool minimums; do not double-count service rates |
| 8. Deliver organism/question/budget recommender with justification and detection limit | `app.py`, `advisor.py`, `web/index.html`, `config/scenarios.json` | Implemented planning advisor. Detection output is analytical read-support probability unless caller truth validation exists | Test feasible and infeasible budgets, custom genome size, rare alleles, repeat goals, insufficient RNA replicates and unavailable empirical evidence |

## Common 60-point rubric

| Criterion | Weight | Reviewable evidence | Current assessment | Residual evidence needed for the strongest descriptor |
|---|---:|---|---|---|
| Biological framing | 5 | Scientific sources; report; recommendation rationale | Implemented and tied to a measured bacterial breadth question | Both students explain why breadth is distinct from variant accuracy, RNA power and assembly quality |
| Databases and acquisition | 8 | Scripted retrieval; measured near-million cohort; manifest; ENA/NCBI joins and conflicts; publication/reference links | Measured, with capped-cohort/full-archive limitation | Complete representative or partitioned archive aggregation would strengthen prevalence claims. ENA/SRA mirrors must not be counted as independent runs |
| Platforms and throughput | 5 | Versioned envelopes, read distributions, costs and references | Implemented planning component | Connect instrument chemistry and throughput with attainable biology; separate historical measurements and current advertised capabilities |
| NGS handling and QC | 7 | Raw/post-audit QC, paired retention 384,691/400,000, adapter/tail/length/quality diagnostics | Measured | Justify thresholds and sampling-limited duplicate diagnostics. A before/after mapping claim needs realignment of trimmed sequences |
| TGS handling and comparison | 7 | Real ONT QC and 50,175/60,000 retention; length distribution, primary gapped alignments and downsampling | Measured; structural assembly/SV advantage remains literature/model-supported | Explain chemistry/library differences and observed edit proxy without an unsupported platform accuracy ranking |
| Alignment core | 8 | Minimap2 2.22 genuine SAM, `sr`/`map-ont`, CIGAR-aware coverage and metrics; separate educational fallback | Measured standard alignment | Explain seeds/chaining/extension, presets, flags and metrics; cross-runtime equality and placement accuracy have not been benchmarked |
| Engineering and scalability | 6 | Streaming acquisition/QC/metadata, deterministic molecule hashes, NumPy difference arrays, measured resource files | Measured bacterial pilot and near-million metadata aggregation | Distinguish transfer, alignment, analysis and aggregation resource scopes. Stored mappings still grow with retained read count |
| Results and interpretation | 6 | Actual curves/QC/cost scenarios, seed ranges, plateau rule, negative finding and limitations | Measured | Present uncertainties, ordered acquisition bias, different BioSamples, price mismatch and unvalidated truth endpoints consistently |
| Reproducibility and code quality | 4 | Small real fixture, replay path, tests, `workflow/Snakefile`, environment/container and README | Implemented; the prepared full native workflow is a different route from the executed pilot | Preserve the checked one-command fixture path and exact environment/tool versions. Full native route requires its own execution before claiming it ran; genuine commits remain required |
| Report | 2 | English report with measured figures, accession table, limitations and AI appendix | Delivered as a separate main-bundle artifact | The report is 10 pages and was rendered and checked. The report specifies task allocation; completed contributions must be supported by the students' actual work |
| Defence and teamwork | 2 | Contribution paragraph, commit history and rehearsed presentation | Required from team | Both named partners make and can explain technical contributions; each can answer questions about every component |

## Honest result language

Use these distinctions consistently in the README, report, interface and presentation:

| If the evidence is… | Say… | Do not claim… |
|---|---|---|
| Ordered streamed portion of an archived read file | “Real acquired prefix of run X; random subsets within this pool” | “Random representative sample of the entire run” |
| Custom exact-anchor interval coverage | “Conservative anchor-supported breadth proxy” | “CIGAR-accurate genome coverage”, “variant sensitivity” or “alignment identity” |
| Independent mapping against a fixed reference | “Mapping fraction remained approximately stable; breadth decreased” when measured | A fabricated decline because the assignment assumes one |
| Highest depth still improving | “Saturation was not reached within the tested pool” | “The largest tested depth is optimal” |
| Poisson/Binomial support calculation | “Probability of at least k supporting observations under stated assumptions” | “Validated limit of detection” |
| Vendor-derived input reproduced by the model | “Arithmetic consistency check” | “Independent validation” |
| Facility sequencing price comparison | “Matching-scope sequencing-component benchmark” | “Validation of the total wet-lab/cloud cost” |
| One organism and historical chemistry | “Bacterial pilot; current-platform recommendations remain conditional” | “Universally validated recommender” |
| An unexecuted extended workflow | “Prepared extension” | “Completed analysis” |

## Submission checks that the team must close

1. Names are filled as Zhaulybaeva Gaukhar and Kurmangazykyzy Dilyara. Confirm actual technical contributions and the course submission date.
2. Check every stated accession, retrieval date and file checksum against the final manifest.
3. The included real small fixture verifies parsing, alignment and the analysis workflow; it does not reproduce the headline deep-pool plateau. Run `python run_experiment.py --pilot` to retrieve/analyse the larger bounded inputs for the headline experiment. Keep the native fastp/SAMtools extension identified separately.
4. Keep large downloaded sequence files out of Git. Include retrieval scripts, small fixture, provenance and checked result summaries instead.
5. Standard minimap2 alignment and deep bounded pools have been executed. Do not describe the separate native fastp/SAMtools workflow or full-archive metadata aggregation as executed unless their own outputs exist. State ordered-prefix, cap and endpoint limits.
6. Confirm report length 8–12 pages, accessible figures, limitations, contribution paragraph and AI disclosure. Retain checked scripts, environment and workflow files.
7. Both partners must make real technical contributions and authentic commits. Do not create backdated commits, impersonate the second partner, or invent a work history.
8. Publish/share the repository with the teaching team, attach the report, and submit before the course deadline. Confirm the actual calendar deadline with the instructor/course page; “week 5” is not a calendar date.
9. Rehearse the ten-minute presentation and five-minute questions with both partners. Explain sampling, flat mapping curves, true versus proxy coverage, platform chemistry, cost scope and the distinction between support probability and detection accuracy.

The remaining items include authentic teamwork, human verification/submission and any broader scientific claims as well as GitHub publication. The core measured pilot is present. A polished repository and presentation cannot replace genuine technical contributions or understanding.
