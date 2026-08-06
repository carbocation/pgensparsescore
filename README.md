# pgensparsescore

`pgensparsescore` scores many polygenic scores in one pass over the relevant
variants of a PGEN file.  Its weight index is variant-major: a decoded dosage
vector is applied only to the scores which contain that variant.

This is an early correctness-first implementation. It supports biallelic
variants, ordinary and conditional-rANS PGEN input, PGEN hardcalls and dosages,
fixed-frequency or cohort-frequency mean imputation, and weights whose effect
allele is either REF or ALT. The allele order in a weight row does not need to
match the PVAR allele order. PGEN dosages are used on the stored 0--2 scale;
chromosome- and sex-specific ploidy transformations are not yet performed.

## Build

The executable links to the PGEN reader libraries in a pinned revision of
[carbocation/plink-ng](https://github.com/carbocation/plink-ng), or an existing
checkout can be supplied for offline builds.

```sh
cmake -S . -B build -DPLINK_NG_SOURCE_DIR=/path/to/plink-ng
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Omit `PLINK_NG_SOURCE_DIR` to let CMake fetch the pinned revision.
Development packages for zlib and zstd are also required.

The same command reads ordinary and conditional-rANS PGENs; storage mode is
detected from the PGEN header. The JSON run metadata records the mode used for
each input file.

## Large score collections

A large collection should be built as one variant index and several score
fragments. The variant index assigns a stable ordinal to each selected
pseudobiallelic variant and recognizes both the score-file ID and the cohort
PVAR ID. This example uses generic column names:

```text
SOURCE_ID	TARGET_ID	REF	ALT
score-id-1	pvar-id-1	A	G
score-id-2	pvar-id-2	C	T
```

```sh
pgensparsescore build-variant-index \
  --variant-list selected-variants.tsv.gz \
  --block-size 100000 \
  --out selected.index.bin
```

Other input column names can be named explicitly with
`--source-id-column`, `--target-id-column`, `--ref-column`, and
`--alt-column`. `--strip-target-id-prefix chr` can, for example, make a
`chr1:...` target alias match a `1:...` PVAR ID without writing a translated
copy of the variant list. The index is a memory-mapped hash table. Building and reading
it does not create a heap-resident set of every variant ID.

Divide the score manifest into groups with similar estimated nonzero-weight
counts. Each group is compiled independently and can run on a different
machine:

```sh
pgensparsescore compile-fragment \
  --manifest fragment_00000.tsv \
  --variant-index selected.index.bin \
  --temp-dir work/fragment_00000 \
  --out fragment_00000.bin
```

A fragment manifest uses the usual `SCORE_ID`, `COLUMN_NAME`, and `PATH`
columns. `SCORE_ID` is the durable identity. Weights store fragment-local
score indices, so an existing fragment remains usable if a later release
adds scores or changes the output column order. The compiler streams weight
rows into variant blocks and sorts one block at a time; it does not retain the
whole fragment in memory.

The fragment list contains one path per compiled fragment:

```text
FRAGMENT
fragments/fragment_00000.bin
fragments/fragment_00001.bin
```

The optional score schema fixes the requested row order and column names:

```text
SCORE_ID	COLUMN_NAME
source:score_a	source__score_a
source:score_b	source__score_b
```

Score all fragments together:

```sh
pgensparsescore \
  --pfile-list cohort.pfiles.tsv \
  --variant-index selected.index.bin \
  --fragment-list fragments.tsv \
  --score-schema score_schema.tsv \
  --read-freq reference.acount.zst \
  --missing-freq omit \
  --out results/cohort
```

The scorer opens every fragment and merges the current variant block from
each one. A PGEN variant is decoded once and its dosage is applied to all
matching weights, regardless of the number of fragments. Memory is bounded by
the index, the current variant's edges, cohort-location and frequency arrays,
and the file-backed score matrix. Fragmentation therefore shortens and
parallelizes the build without multiplying genotype reads. During scoring,
the binary divides each dosage's samples among a persistent worker pool using
all physical cores visible to the process. `--threads N` overrides that
default. Additional threads do not read the same genotype again.

## Monolithic compiled catalogs

The manifest form is convenient for small runs, but it requires every score
file to be parsed each time a cohort is scored. A compiled catalog performs
that parsing once and stores the nonzero weights grouped by source variant.
It is independent of a particular PGEN and can be reused for 1KG+HGDP, AoU,
UKB, and other cohorts.

This format is retained for smaller collections and compatibility. The
variant-index and fragment format above is the intended path when all weights
cannot be built comfortably in one process.

```sh
pgensparsescore compile \
  --manifest scores.tsv \
  --variant-map selected-variants.tsv.gz \
  --out selected.catalog.bin
```

In compile mode, the `SOURCE_ID` keys in `--variant-map` select the variants
included in the catalog; its `TARGET_ID` values are not stored. Weight-zero
rows are recorded in QC counts and omitted from the catalog. Without a variant
map, every nonzero weight in the manifest is included.

Repeated rows for the same variant within one score are retained, so their
contributions add just as they do in the source score. Rows after the first are
also counted as duplicates in the score metadata.

Score a cohort by supplying the compiled catalog in place of the manifest:

```sh
pgensparsescore \
  --pfile-list cohort.pfiles.tsv \
  --compiled-catalog selected.catalog.bin \
  --variant-map cohort.variant-map.tsv.gz \
  --read-freq reference.acount.zst \
  --missing-freq omit \
  --out results/cohort
```

The binary format stores score metadata, source variant IDs, unordered allele
pairs, and effect weights. REF/ALT orientation is still determined and checked
against each target PVAR when the catalog is used. The companion `.json` file
records its dimensions.

## Inputs

The compact score manifest is a tab-separated file with `SCORE` and `PATH`
columns. Relative paths are resolved relative to the manifest.

```text
SCORE\tPATH
score_a\tweights/score_a.tsv.gz
score_b\tweights/score_b.tsv.gz
```

Compiler manifests may instead provide `SCORE_ID`, `COLUMN_NAME`,
`DISPLAY_NAME`, and `PATH`. `COLUMN_NAME` is then used for the score-table
column. This keeps the source identifier and descriptive label in the manifest
while using an analysis-safe name in Parquet and R. Older `SCORE` and
`SCORE_ID` manifests remain accepted.

Each weight file is a tab-separated plain-text, gzip-compressed, or
zstd-compressed file with these columns:

```text
SNP\tEFFECT_ALLELE\tOTHER_ALLELE\tEFFECT_ALLELE_WEIGHT
```

`SNP` is matched exactly to the PVAR `ID`.  The effect and other alleles must
be the PVAR REF/ALT pair, but they may appear in either order.  For example,
both of these are valid for a PVAR record with `REF=C` and `ALT=T`:

```text
EFFECT_ALLELE=T  OTHER_ALLELE=C
EFFECT_ALLELE=C  OTHER_ALLELE=T
```

The scorer does not perform strand-complement guessing.  A matched variant
whose two alleles do not exactly equal the PVAR REF/ALT pair is an error.
Multiallelic PVAR records are rejected when referenced by a score.

### Variant-ID mapping

When the weight files and target PVAR use different IDs for the same exact
REF/ALT variant, supply a two-column variant map:

```text
SOURCE_ID\tTARGET_ID
DRAGEN:chr1:100:A:G\trs123
```

```sh
pgensparsescore \
  --pfile-list reference.pfiles.tsv \
  --manifest scores.tsv \
  --variant-map reference.variant-map.tsv.zst \
  --read-freq reference.acount.zst \
  --error-on-missing-freq \
  --out results/reference
```

The map may be plain text, gzip, or zstd. `SOURCE_ID` is the ID in each score
row and `TARGET_ID` is the PVAR ID. Both columns must be unique. A score row
whose source ID is absent from the map is counted as a missing variant; a
mapped row is still required to have the exact PVAR REF/ALT allele pair. This
keeps ID translation separate from allele interpretation and avoids creating
a translated copy of every score file.

When a map is supplied, the scorer scans each PVAR but retains metadata only
for its `TARGET_ID` values. A compiled catalog also limits retained PVAR rows
to variants present in that catalog, even when no ID map is needed. PGEN reads
still use the original PVAR row numbers. This permits scoring a small catalog
from a much larger PGEN without first writing a subset PGEN.

## Run

```sh
pgensparsescore \
  --pgen cohort.pgen \
  --pvar cohort.pvar \
  --psam cohort.psam \
  --manifest scores.tsv \
  --out results/catalog
```

To score variant-sharded PGENs, such as one file per chromosome, provide a
tab-separated PGEN list instead:

```text
PGEN\tPVAR\tPSAM
cohort.chr1.pgen\tcohort.chr1.pvar.zst\tcohort.chr1.psam
cohort.chr2.pgen\tcohort.chr2.pvar.zst\tcohort.chr2.psam
```

```sh
pgensparsescore \
  --pfile-list cohort.pfiles.tsv \
  --manifest scores.tsv \
  --out results/catalog
```

Relative paths in the PGEN list are resolved relative to that list.  Every
PGEN must represent a distinct set of variants and have exactly the same
samples, FIDs when present, and sample order.  The score manifest is compiled
once against the union of the PVARs, then partitioned by PGEN.  PGENs are
decoded sequentially into one file-backed accumulator and the output table is
written only once.

### Fixed-frequency imputation

Use `--read-freq` when the same reference frequency must impute missing
genotypes in more than one target dataset:

```sh
pgensparsescore \
  --pfile-list cohort.pfiles.tsv \
  --manifest scores.tsv \
  --read-freq reference.acount.zst \
  --error-on-missing-freq \
  --out results/catalog
```

The frequency file can be plain text, gzip, or zstd.  It must have `ID`,
`REF`, and `ALT` (or `ALT1`), plus one of the following representations:

- `ALT_DOSAGE_MEAN`;
- `ALT_FREQS` or `ALT1_FREQ`;
- `ALT_CTS` or `ALT1_CT`, together with `OBS_CT`; or
- `REF_CT`, together with `OBS_CT`.

This includes biallelic PLINK 2 `.acount[.zst]` and `.afreq[.zst]` files.  An
allele frequency `p` is converted to expected ALT dosage `2p`.  For every
scored variant found in the frequency table, the frequency REF and ALT must
exactly match the PVAR REF and ALT; disagreement is always an error.

`--missing-freq` controls a matched PGEN variant which has no frequency row:

- `cohort` uses its nonmissing mean in the PGEN being scored. This is the
  default and preserves the earlier behavior.
- `error` stops the run. The older `--error-on-missing-freq` spelling remains
  accepted.
- `omit` removes the variant from every score, including any REF-effect
  intercept. This is useful when scores must be restricted to variants
  represented by a normalization reference.

With `--read-freq`, only frequency rows relevant to matched score variants are
retained in memory. The input file is still checked as a stream.

The primary output, `catalog.scores.tsv.gz`, is a gzip-compressed wide table.
Each row is a sample, followed by one column named for every manifest `SCORE`.
`IID` is required in the PSAM.  If the PSAM also has `FID`, it is preserved;
the scorer does not invent one when it is absent:

```text
IID\tscore_a\tscore_b
sample1\t0.125\t-1.75
sample2\t0.5\t0.25
```

Additional outputs are:

- `catalog.score-metadata.tsv`: score order, weight inclusion, matching,
  repeated variants, frequency omission, and allele-orientation counts.
- `catalog.json`: output dimensions and run-level counters.

The JSON also records `variant_mapping_rows` (zero when `--variant-map` is not
used), `pvar_variants_loaded`, and the storage mode of each input PGEN.

## Progress reporting

Both catalog compilation and scoring can write a structured progress log:

```text
--progress-jsonl run.progress.jsonl --progress-interval-seconds 30
```

The JSON Lines file is flushed after every event. It records the current phase,
elapsed time, current and peak process memory when the operating system exposes
them, and phase-specific counters. Compilation reports score files, input and
retained weights, unique variants, sorting, and serialized bytes. Scoring
reports catalog loading and matching, PVAR and frequency loading, sparse and
dense decodes, weight edges assigned to each decode path, estimated sparse and
dense score updates, imputed values, and output rows. Phase-boundary events are
always written; long phases also write events at the requested interval.

Variant-index construction reports both passes through the list, hash-table
size, and output bytes. Fragment construction reports files and weights read,
excluded and repeated rows, blocks serialized, and output bytes. Fragment
scoring reports fragment bytes opened, variant groups merged, and total PGEN
decodes; the decode count is independent of the fragment count. It also reports
the resolved thread count and how many score updates used the parallel kernel.

The scorer uses a file-backed score-major matrix while applying variants,
since that is the efficient update layout.  It transposes that working matrix
in bounded-memory blocks when writing the sample-major table, then removes the
working file.  Score-major layout is therefore an implementation detail, not
part of the output contract.

## Scoring semantics

The input allele order determines which of two equivalent calculations is
used; it is not an input restriction:

| Weight row | Contribution |
| --- | --- |
| effect = PVAR ALT, other = PVAR REF | `w * ALT_DOSAGE` |
| effect = PVAR REF, other = PVAR ALT | `w * (2 - ALT_DOSAGE)` |

Internally, the second expression is rewritten as
`2 * w - w * ALT_DOSAGE`.  This permits the scoring kernel to operate entirely
on ALT dosages without changing the meaning of a REF-effect weight.  Missing
ALT dosages are mean-imputed from the external expected ALT dosage when
`--read-freq` supplies the variant.  Otherwise they are imputed from nonmissing
samples at that variant.  An all-missing scored variant is permitted with an
external frequency and is an error without one.

When `pgenlib` can return a common dosage plus a short difference list, the
common contribution is accumulated once per score and only the differing
samples are updated.  Otherwise, the dense dosage vector is decoded and
applied.  The score matrix is file-backed, so RAM does not grow as
`sample_count * score_count`.

## License

The application is GPL-3.0-only.
