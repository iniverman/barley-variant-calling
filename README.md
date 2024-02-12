## Introduction


In this repository we perform a variant calling study using two reads of a barley sample. For that we used freebayes to perform the variant study. 
This is part of a project of practices done in EEAD-CSIC 

Variant calling allow us to detect and study single nucleotide polymorphisms (SNPs) and insertions and deletions (indels) ,among other, from a sequence previously get from next generation sequencing (NGS). To do this we will use freebayes, a bayesian genetic variant detector that is haplotype based, that means that it calls variant based on the literal sequenced reads aligned to a particular target, not their precise alignment.


## Materials & Methodologies


To perform this study we will use Snakemake from anaconda and the freebayes and picard tools. 

We will use a pired end couple of samples of bareley (A_1_20_1,A_1_20_2), in this case both of them are in format .fastq.gz. Also we use as reference the barely genome (GCA_904849725.1_MorexV3_pseudomolecules.chrnames.fna).


The steps to do out variant calling are:

  1- Preparation of our data.

  2- Quality Control of our samples.

  3- Mapping against the reference genome.
  
  4- Mark duplicates using Picard.

  5- Perform the variant calling using freebayes. 
  

  ### 1-Preparation of our data.

 
It is possible that the data is compressed (as in our case), so for that we make a rule to descompress our archives. We will use:
  ```
  snakemake data/samples/A_1_20_{1,2}.fastq -c4
  ```
Then we will get the fastq for each of our samples. Fastq is a text format to store biological sequences,obtained from the sequenciation, and their quality scores. 


### 2- Quality Control of our samples.


To get the quality of each of the fastq you have to run this code:
```
data/samples/A_1_20_{1,2}_fastqc.{html,zip}
```
Quality control is important because it give us specific information about the quality of the fastq to study, and can help us deciding if we have to do a preprocessing to our sample, like removing the adapters used in the sequencing for example.


### 3- Mapping against the reference genome.

The first step is to index the reference genome. To do that you can run this code: 
```
bwa index -a bwtsw data/GCA_904849725.1_MorexV3_pseudomolecules.chrnames.fna
```




