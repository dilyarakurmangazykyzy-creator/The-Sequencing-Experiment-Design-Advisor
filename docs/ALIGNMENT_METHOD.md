# Alignment and coverage interpretation

## What the pilot measures

The primary experiment uses real archived Illumina mates and ONT reads against
the same bacterial reference. Minimap2 produces genuine SAM base-level
alignments. The primary pilot records the executable/version in its run
provenance; the execution orchestrator uses the WebAssembly build of minimap2
2.22. This is the actual aligner compiled to another runtime, rather than a
custom mapper given its name. The separately supplied educational 13-mer mapper
is a fallback and has a different, explicitly labelled conservative endpoint.

Illumina uses `-ax sr` with both mate files; ONT uses `-ax map-ont`. Both presets
are documented by the [minimap2 project](https://github.com/lh3/minimap2).
Minimap2 finds representative k-mer seeds, chains consistent seeds to propose
placements, and extends/verifies candidate alignments with dynamic programming.
Presets change this search and scoring for short reads or indel-prone long reads.
The current project also offers newer high-quality ONT presets, but the actual
2.22 execution is intentionally reported with its available `map-ont` preset;
current documentation is not evidence that a newer preset was run.

The biological endpoint is **reference breadth and depth**, not variant recall,
expression accuracy or assembly completeness. A high breadth value establishes
that positions receive accepted alignment evidence; it does not establish that
all true variants can be called or repeats can be assembled.

## Why one alignment pass supports exact nested subset coverage

The pipeline normalizes the molecule identifier and hashes it with BLAKE2b and
an explicit seed. The same identifier selects both Illumina mates. Fractions
are nested: every selected molecule at 10% is selected at 20% with that seed.
Seeds 17, 29 and 43 create three alternative molecule subsets of the acquired
data. The randomization acts within the acquired FASTQ prefix, not across the
unavailable complete archived run.

Let `I(u)` be the accepted SAM reference intervals of molecule `u` and `h(u,s)`
its hash fraction at seed `s`. At subset fraction `f`, depth at position `p` is:

```text
d(p,f,s) = sum over accepted molecules u with h(u,s)<f
           of the number of SAM intervals in I(u) containing p
```

Computing this sum from stored intervals is exactly equivalent to taking the
same molecule subset of the produced SAM and computing its coverage under the
same flag/MAPQ/CIGAR policy. It avoids re-running the expensive alignment 27
times. The reference, input sequences, preset and paired relationship are fixed.
Ordinary reference mapping evaluates each read or mate pair rather than using
the sequencing depth of other molecules to improve its placement.

This exact claim concerns a subset of **the produced SAM**. A separate empirical
comparison of fresh native re-alignments and the WebAssembly execution has not
been performed. Threading, implementation versions, batching, changed trimming,
reference choice or scoring parameters can change alignment output. Therefore
the pilot does not claim a measured cross-runtime equivalence benchmark.

The optional `workflow/Snakefile` runs fastp, native minimap2 and samtools for
each subset, using the pinned environment and per-alignment resource logs. This
full route must be identified as unexecuted until its actual output exists.
Its trimming is different from the submitted raw-alignment audit; its results
must not be substituted silently for the pilot's result.

## QC selection and the raw SAM

The pilot performs a streaming raw-read QC audit and tests platform-specific
eligibility. Illumina trimming is an audit of the common 12-base adapter prefix
and low-Q20 tails, followed by minimum-length/mean-quality/N filters. If either
mate fails, both are rejected. ONT uses minimum length500 and mean Phred7, with
no Illumina adapter or short-read tail assumption. Thresholds and rejection
counts are recorded rather than inferred from an ideal platform specification.

The external SAM corresponds to the acquired **raw** FASTQ. Accepted molecule
IDs are selected from that SAM; internally audited trimming does not re-align
their sequences. SAM intervals can therefore still include a trimmed tail if
the aligner chose to align it, or omit it if soft-clipped. Nominal depth uses the
original lengths of retained raw molecules so its numerator corresponds to the
actual external alignment input. Post-filter QC lengths describe the audit and
must not be described as the sequences supplied to the external aligner.

This scope is an accepted-molecule analysis of raw alignments. A claim that
adapter trimming improved mapping would require an actual before/after
alignment experiment, which this pilot does not provide.

## Coordinate, CIGAR and filtering rules

SAM positions are one-based; the coverage code subtracts1 and stores zero-based
half-open intervals `[start,end)`. CIGAR `M`, `=` and `X` contribute reference
positions. Insertions consume query bases without reference depth. Deletions
and skipped intervals advance the reference coordinate but provide no covered
bases. Soft/hard clipping does not contribute to reference depth. Adjacent match
intervals across an insertion are merged without changing any covered position.
Contig bounds are checked before coverage is calculated.

Reference breadth uses mapped primary alignments with MAPQ≥20 and excludes
secondary, supplementary and preflagged QCFAIL/DUP records. MAPQ255 means
unavailable and is excluded. Duplicate marking has not been run, so unflagged
PCR/biological duplicates remain. Both mates contribute independently where they
overlap; depth is aligned-read depth, not independent-fragment depth.

The mapped-read numerator counts all mapped primary reads among retained
molecules, including low-MAPQ reads and preflagged QCFAIL/DUP records. Its
denominator counts all their primary reads. Consequently a high mapping rate
does not guarantee high MAPQ-filtered breadth. MAPQ≥20 is an operational
confidence filter, not proof of unique correct placement. Minimap2's placement
score is heuristic and has not been independently probability-calibrated here.
The [official option reference](https://lh3.github.io/minimap2/minimap2.html)
explains presets, scoring and SAM output options.

The reported `mean_alignment_edit_distance_per_aligned_query_base` averages
`NM/(M+I+=+X)` across records with an NM tag. NM combines substitutions and
indels, including both sequencing errors and biological reference differences.
Deletions enter NM but not that query-base denominator. This ratio is an edit
distance proxy, **not** sequencing-error truth or exact sequence identity.

## Depth-dependent results and limitations

The pipeline reports nominal bases/genome length, mapped mean depth, breadth at
≥1×/5×/10×, per-reference-window GC and dropout, and mapping rate. The three-seed
minimum/maximum range describes sampling variability within the acquired data;
it is not a confidence interval over biological samples or archived studies.
The Poisson `1-exp(-depth)` curve is an ideal uniform independent-sampling
baseline. Real GC bias, overlapping mates, duplicates, soft clipping and low-MAPQ
regions can cause deviations.

The descriptive plateau requires at least95% breadth≥1× and two adjacent gains
under one percentage point. It is specific to the chosen grid, reference,
filters and available depth. Failure to reach it is recorded as `null`. Even a
detected plateau does not establish sufficient depth for variant calling,
expression or assembly.

With fixed alignments, removing unrelated molecules does not causally degrade
the placement quality of a surviving read. Small subset mapping-rate changes
reflect which reads were selected; genomic breadth and evidence per position
are the quantities that fall with depth. A flat mapping-rate curve is a useful
negative result, rather than a failed experiment.

The two runs are separate experiments and BioSamples; sequencing chemistry,
library preparation including whole-genome amplification, archived read order
and study conditions confound their comparison. Prefix sampling, one bacterial
reference, missing truth calls and missing independent biological replication
prevent a general platform accuracy ranking or transfer to human/RNA/assembly
requirements. Recommendations for those applications remain model assumptions
with explicit validation gaps.

## Engineering and measured scale

FASTQ QC and metadata aggregation stream input records. Retained per-molecule
SAM intervals still occupy memory, and coverage uses NumPy difference arrays or
a standard-library sparse sweep. This is a bacterial pilot, not an assertion
that arbitrary human multi-terabyte SAM files fit in RAM. The full workflow uses
sorted BAM and streamed samtools depth instead.

`summary.json` records pilot processing time and peak memory; alignment runtime
belongs to the separate actual aligner logs. `metadata_summary.json` records
actual row count, parse/aggregate runtime, rows per second and peak process
memory. `metadata_scale.json` profiles ordered1000/10000/100000-row prefixes in
fresh processes and includes the measured full acquisition. Ordered prefixes
are valid computational loads, not probability samples for design frequencies.
Any linear10million-row estimate is labelled a projection, not a completed run.
`qc_alignment_comparison.tsv` places both platforms side by side with physical
raw-read counts, actual read lengths, raw and audit-quality statistics, accepted
identifier counts, mapping/MAPQ/edit-distance proxies and filtered breadth.
Every row includes its unit and scope so read pairs, raw-SAM input and internally
audited trimmed sequences are not conflated.

Archive `base_count/read_count` ratios are **bases per returned archive count
unit**. In the selected paired SRR1030394 run, read_count agrees with pair/spot
count, so the ratio is roughly twice the physical mate length. Other deposited
layouts need not share that convention; all PAIRED rows are not blindly divided
by2. Physical mate lengths come from FASTQ QC. Fields absent from the requested
TSV schema are separated from returned columns containing blanks; schema
absence is not evidence of bad archive metadata. The acquisition follows the
[official ENA Portal API documentation](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/advanced-search.html).
