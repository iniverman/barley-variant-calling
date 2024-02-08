## Introduction

In this repository we perform a variant calling study using two reads of a barley sample. For that we used freebayes to perform the variant study. 
This is part of a project of practices done in EEAD-CSIC 

Variant calling allow us to detect and study single nucleotide polymorphisms (SNPs) and insertions and deletions (indels) ,among other, from a sequence previously get from next generation sequencing (NGS). To do this we will use freebayes, a bayesian genetic variant detector that is haplotype based, that means that it calls variant based on the literal sequenced reads aligned to a particular target, not their precise alignment.

## Materials & Methodologies
To perform this study we will use Snakemake from anaconda and the freebayes and picard tools.

The steps to do out variant calling are:

  1- Preprocessing and preparation of our data.

  2- Quality Control of our samples.

  3- Map the reads against the reference genome.
  
  4- Mark duplicates using Picard.

  5- Perform the variant calling using freebayes. 

First we will unzip our sequences as they are compressed, for that we will use:
   ```
  snakemake data/samples/A_1_20_{1,2}.fastq -c4
  ```



