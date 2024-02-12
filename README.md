## Introduction


In this repository we perform a variant calling study using two reads of a barley sample. For that we used freebayes to perform the variant study. 
This is part of a project of practices done in EEAD-CSIC 

Variant calling allow us to detect and study single nucleotide polymorphisms (SNPs) and insertions and deletions (indels) ,among other, from a sequence previously get from next generation sequencing (NGS). To do this we will use freebayes, a bayesian genetic variant detector that is haplotype based, that means that it calls variant based on the literal sequenced reads aligned to a particular target, not their precise alignment.


## Materials & Methodologies


To perform this study we will use Snakemake from anaconda, samtools, bwa and the freebayes and picard tools.

We will use a pired end couple of samples of bareley (A_1_20_1,A_1_20_2), in this case both of them are in format .fastq.gz. Also we use as reference the barely genome (GCA_904849725.1_MorexV3_pseudomolecules.chrnames.fna).


The steps to do out variant calling are:

  1- Preparation of our data.

  2- Quality Control of our samples.

  3- Mapping against the reference genome.

  4- Marking duplicates and treating the readgroups.

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
snakemake data/samples/A_1_20_{1,2}_fastqc.{html,zip} -c4
```
Quality control is important because it give us specific information about the quality of the fastq to study, and can help us deciding if we have to do a preprocessing to our sample, like removing the adapters used in the sequencing for example.


### 3- Mapping against the reference genome.

The first step is to index the reference genome. To do that you can run this code: 
```
bwa index -a bwtsw data/GCA_904849725.1_MorexV3_pseudomolecules.chrnames.fna
```
This is a step that can take a lot of time (from hours to days depending of the genome and the computational capacity) so be patient. You can estimate the time is going to take because the terminal will show you so often how many characters have it study so far.

Then we can proceed with the mapping of our samples against the genome reference:
```
snakemake results/mapped/A_1_20.bam -c4
```
The BAM file is the compressed binary version of a SAM file,that is the format to save the mapping of the sequences.
Then we sort the BAM file wih:
```
snakemake results/sorted/A_1_20_sorted.bam -c4
```

#### 3.5- Calculate and visualize the depth of the mapping.

It is important to check how the map process has gone. To do that we check the depth of our mapped file with: 
```
snakemake results/depth/A_1_20_depth.csv -c4
```
With this we can check the depth for chromosome and position in the chromosome, we can even check the mean of the depth for each one using: 
```
snakemake results/depth/A_1_20_mean.csv -c4
```
But to see it in a more easy and understandable way, we can print it graphically and make it more visual with this:

```
snakemake results/depth/plots/A_1_20.svg -c4

```
Using this we will get the plots for the depth for each chromosome. 

### 4- Marking duplicates and treating the readgroups.

Now we are going to mark the duplicated reads caused during the sequenciation step, for example during the PCR. To do this we use:

```
snakemake results/markdup/A_1_20_nodup.{bam,txt} -c4
```







