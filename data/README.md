# Data

Inputs for `pathogenomics run`. Everything here is either public data (NCBI, wwPDB,
AlphaFold DB) or output I generated from public data in Galaxy (usegalaxy.org) during
BIOL 4640 in Fall 2025. Files were renamed to stable, descriptive names and large text
files were gzip-compressed (`gzip -9 -n`); contents are otherwise unchanged. Course
handouts, slides and the syllabus are the instructor's property and are not included.

## Provenance

| Path | What it is | Original name / source |
| --- | --- | --- |
| `qc/reads/BR-Tn1_S4_L001_R1_001.fastq.gz` | 5,190 single-end Illumina reads (78 bp) | Activity 4 input |
| `qc/fastqc/BR-Tn1_R1.fastqc_data.txt.gz` | FastQC 0.12.1 on the reads above | Galaxy "FastQC on data 45" |
| `qc/fastqc/EN59_R1.fastqc_data.txt.gz`, `EN59_R2…` | FastQC on EN59_S48 R1 and R2 | Galaxy "FastQC on data 39/40" |
| `qc/fastqc/SRR29263899_R1.fastqc_data.txt.gz` | FastQC on SRR29263899 forward reads | Galaxy "FastQC on data 31" |
| `qc/fastqc/SRR10828193_R1.fastqc_data.txt.gz` | FastQC on SRR10828193 forward reads | Galaxy "FastQC on data 21" |
| `qc/fastqc/activity3_reverse.fastqc_data.txt.gz` | FastQC on the Activity 3 reverse reads | Galaxy "FastQC on data 8" |
| `qc/cutadapt/*.cutadapt.txt` | Cutadapt 5.1 reports (parameters on line 2 of each) | Galaxy datasets 55, 58, 65, 72 |
| `assembly/pseudomonas_aeruginosa/SRR1213109{3,4,6,7}.spades_contigs.fasta.gz` | SPAdes contigs, 4 *P. aeruginosa* isolates | Final exam Q42 ("SPAdes on dataset … Contigs") |
| `assembly/pseudomonas_aeruginosa/reference_genomes.tsv` | sizes of 5 complete *P. aeruginosa* RefSeq genomes | NCBI Datasets `data_summary.tsv` (midterm download) |
| `assembly/staphylococcus_aureus/activity6_genome.fasta.gz` | single-sequence *S. aureus* genome annotated in Activity 6 | `Activity6_Section2_genome.fasta` |
| `assembly/outbreak_sample/outbreak_sample.shovill_skesa.fasta.gz` | Shovill 1.1.0 / SKESA assembly (42 contigs, 4.97 Mb, 52% GC – not *S. aureus*) | `Outbreak_sample.fasta` |
| `annotation/staphylococcus_aureus/refseq/GCF_*.gff.gz` | RefSeq annotation of 6 complete *S. aureus* genomes | NCBI Datasets (midterm) |
| `annotation/staphylococcus_aureus/refseq/strains.tsv` | accession → strain, size, gene count | NCBI Datasets `data_summary.tsv` |
| `annotation/staphylococcus_aureus/bvbrc/activity6_genome.bvbrc.gff.gz` | BV-BRC (RASTtk) annotation of the Activity 6 genome | `Staphylococcus aureus Activity6.gff` |
| `annotation/staphylococcus_aureus/bvbrc/activity6_genome.genome_quality_details.txt` | BV-BRC completeness/contamination | `genome_quality_details.txt` |
| `pangenome/roary_gene_presence_absence.csv.gz` | Roary gene presence/absence: Prokka-annotated EN93, EN94, EN96, EN97 + PA14 | Galaxy dataset 168 (final) |
| `rnaseq/mabs/counts/SRX2122277{1..6}.featureCounts.tsv` | featureCounts (Bowtie2 alignments) vs *M. abscessus* ATCC 19977 | Activity 10 |
| `rnaseq/mabs/samples.tsv` | design: SRX21222771–73 = "condition", 74–76 = "WT R" (wild-type rough) | Galaxy sample labels |
| `rnaseq/mabs/galaxy_deseq2_results.tsv` | DESeq2 results from Galaxy (header added) | `Galaxy19.xlsx` |
| `rnaseq/mabs/galaxy_normalized_counts.tsv` | DESeq2 normalised counts from Galaxy | Galaxy dataset 16 |
| `rnaseq/mabs/GCF_000069185.1.gff.gz` | RefSeq annotation of *M. abscessus* ATCC 19977 | NCBI Datasets (Activity 9) |
| `rnaseq/mtb/counts/SRX368674{9..54}.featureCounts.tsv` | featureCounts vs *M. tuberculosis* H37Rv | Final exam Q43 |
| `rnaseq/mtb/samples.tsv` | design: SRX3686749–51 = control, 52–54 = experimental | final exam design |
| `rnaseq/mtb/galaxy_deseq2_results.tsv` | DESeq2 "experimental vs control" from Galaxy (header added) | Galaxy dataset 89 |
| `rnaseq/mtb/GCF_000195955.2.gff.gz` | RefSeq annotation of *M. tuberculosis* H37Rv | NCBI Datasets (final) |
| `structure/alphafold3_zur_Rv2359/` | AlphaFold Server (AlphaFold 3) job for Mtb Zur/FurB (Rv2359): request, 5 summary confidences, model 0 (mmCIF) and its full data (pLDDT, PAE) | `fold_2025_12_08_18_11*` |
| `structure/afdb/AF-B1MH84-F1-model_v6.pdb` | AlphaFold DB model, *M. abscessus* LuxR-family regulator | alphafold.ebi.ac.uk |
| `structure/afdb/AF-P9WMF9-F1-model_v6.pdb` | AlphaFold DB model, *M. tuberculosis* DevR/DosR | alphafold.ebi.ac.uk |
| `structure/pdb/6vxx.pdb.gz`, `9f9y.pdb.gz` | SARS-CoV-2 spike: Wuhan-Hu-1 closed state; BA.2.87.1 ectodomain | rcsb.org (`pdb6vxx.ent`, `pdb9f9y.ent`) |

Raw reads not redistributed here (too large) can be fetched from the SRA by the accessions
above (e.g. `prefetch SRR12131093`).

## Terms

* **NCBI** (RefSeq annotation, assembly reports, SRA-derived data): public, see
  <https://www.ncbi.nlm.nih.gov/home/about/policies/>.
* **wwPDB** entries 6VXX and 9F9Y: CC0 1.0.
* **AlphaFold DB** models: CC BY 4.0 – Jumper *et al.*, *Nature* 2021; Varadi *et al.*,
  *Nucleic Acids Res.* 2024.
* **AlphaFold Server output** (`structure/alphafold3_zur_Rv2359/`): provided under and subject to
  the [AlphaFold Server Output Terms of Use](https://alphafoldserver.com/output-terms), included as
  `TERMS_OF_USE.md` – non-commercial use only; not to be used to train structure-prediction models.
* Galaxy outputs derived from the data above: my own work, shared under the repository's MIT licence.
