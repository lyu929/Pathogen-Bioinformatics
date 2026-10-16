# Methods

## What was run in Galaxy (course work)

The upstream analyses were performed on <https://usegalaxy.org> as part of the weekly activities,
the take-home midterm and the final exam. Parameters below are the ones recorded in the
output files themselves.

| Step | Tool | Recorded parameters |
| --- | --- | --- |
| Read QC | FastQC 0.12.1 | defaults |
| Trimming | Cutadapt 5.1 | `--quality-cutoff=20 --minimum-length=22 --error-rate=0.1 --overlap=3 --times=1`; `-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` for the Rough_SCFM_3 small-RNA library |
| Assembly | SPAdes (Galaxy); Shovill 1.1.0 with SKESA for the outbreak sample | contig headers record SPAdes k-mer coverage / Shovill read coverage |
| Annotation | Prokka (Galaxy) for the pangenome isolates; BV-BRC/RASTtk for the Activity 6 genome; NCBI RefSeq (PGAP) for the reference genomes | – |
| Pangenome | Roary | Galaxy defaults |
| Alignment / counting | Bowtie2 → featureCounts | gene-level counts on the RefSeq annotation |
| Differential expression | DESeq2 (Galaxy) | `~ condition`; *M. abscessus*: condition vs WT rough; *M. tuberculosis*: experimental vs control |
| Structure | AlphaFold Server (AlphaFold 3), AlphaFold DB v6 models, RCSB PDB | AF3: one protein chain, 1 seed, templates on |

## What this repository computes

* **Quality trimming** – `qc.quality_trim_index` implements the modified BWA algorithm that
  Cutadapt uses for `--quality-cutoff`: scanning from the 3′ end, accumulate `cutoff − q` and
  cut at the position where the running sum peaks; stop when it turns negative. The 5′ cutoff is
  the mirror image. Reads shorter than `--minimum-length` after trimming are discarded.
* **FASTQ statistics** – per-position mean Phred quality (pooled over all reads reaching a
  position, as FastQC does for its base groups), GC% over A/C/G/T, length range.
* **Assembly metrics** – Nx/Lx (and NGx when an expected genome size is given) following the
  QUAST definitions; GC; N density; length-weighted mean coverage parsed from contig headers.
* **Annotation summary** – feature counts, pseudogenes (`gene_biotype=pseudogene` or
  `pseudo=true`), coding density as the union of non-pseudo CDS intervals over genome length,
  share of hypothetical proteins. Markers are matched on protein *product* text with
  curated regular expressions (`annotation.STAPH_MARKERS`).
* **Pangenome** – Roary partitions (core ≥ 99%, soft-core 95–99%, shell 15–95%, cloud < 15% of
  genomes), permutation rarefaction curves (500 random genome orders, seed 1), Jaccard distances on
  accessory clusters.
* **RNA-seq** – median-of-ratios size factors; log2(normalised + 1) PCA on the 500 most variable
  genes; Wald tests with PyDESeq2 (`~condition`, Cook's-distance filtering and independent
  filtering on, α = 0.05). Results are compared with the Galaxy DESeq2 tables on base means,
  log2 fold changes, p-value ranks and the overlap of significant genes (padj < 0.05,
  |log2FC| ≥ 1).
* **Structure** – pLDDT read from the B-factor column (CA atoms); PAE from AlphaFold 3
  `full_data_0.json`; residue correspondence from global BLOSUM62 alignment (gap −10/−0.5, free end
  gaps) of SEQRES sequences when available; optimal superposition by the Kabsch algorithm
  (reflections excluded); rigid bodies by repeatedly extracting the largest pair set that
  superposes within 2 Å, seeded from every contiguous 12-residue fragment.

## Reproducing

```bash
pip install -e ".[dev]"
pathogenomics -v run      # regenerates results/ and reports/
pytest                    # unit + validation tests
```
