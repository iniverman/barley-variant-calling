# Variant Calling 

## Introduction


In this repository we are going to perform a variant calling study using two samples of barley. To do it we used freebayes,GATK and deepvariant to perform the study. 


With the variant calling we want to detect and study single nucleotide polymorphisms (SNPs), and insertions and deletions (indels) ,among other, from a sequence previously got from next generation sequencing (NGS). 


## Materials & Methodologies


To perform this study we will use Snakemake from anaconda, samtools, bwa, bcftools and the deepvariant,freebayes,picard and GATK tools.  
The data samples are in data/samples, and the genome reference in fasta format is ,in the data folder. The reference genome in this repository is much smaller than the real one so I recommend downloading the whole from here https://www.ebi.ac.uk/ena/browser/view/GCA_904849725.1.

In this case we will use a pired end couple of samples of bareley (A_1_20_1.fastq.gz,A_1_20_2.fastq.gz),both of them are in format fastq.gz. Also we use as reference the barely genome GCA_904849725.1_MorexV3_pseudomolecules.chrnames.fna.

Before starting, please clone this repository to use it :

```
git clone https://github.com/iniverman/barley-variant-calling
```

## Installation guide 

First we need to install all the proper programs that we will use in the work, most of them are included with snakemake and the conda enviroment but others must be installed.

We are going to start installing snakemake, to do it you should see the installation guide and tutorials from the official web: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html. 

If you don't want, you can run directly:

```
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```

Then you need to activate the conda enviroment, for that we create a working directory, for example:

```
mkdir snakemake-tutorial
cd snakemake-tutorial
```

And we activate conda:

```
conda activate base
```

Once we have installed the conda enviroment and we have activated it, we can install the mamba command:

```
conda install -n base -c conda-forge mamba
```

And use it to create the enviroment:

```
mamba env create --name snakemake-tutorial --file environment.yaml

```

Finally we activate the enviroment with:

```
conda activate snakemake-tutorial

```

If you want to deactivate it you need to use :

```
conda deactivate

```

Please take into account that the name of this enviroment is for this example and you can change it, also you should check https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html for a proper tutorial on how to use snakemake.


Now we are going to download and move to the working directory the gatk and picard tools:

1 - For gatk go to https://github.com/broadinstitute/gatk/releases and download release: gatk-4.5.0.0.zip (if it's not the version 4.5.0.0 , the program will not     work).
    
2 - For picard go to https://broadinstitute.github.io/picard/. Then ,in that page, in the top right ,click on Latest Jar Release. Finally scroll down in the new       window and download picard.jar.
    
3 - Copy both of the downloads to the folder 'barley variant calling' that have been created by cloning this repository. Both of the programs must be in the same      folder as the snakefile.


Now we are going to install freebayes.

Once you have installed conda,and you have the enviroment active you can use :

```
conda install bioconda::freebayes
```
to directly install freebayes

You can check the documentation for freebayes in https://github.com/freebayes/freebayes




## Variant Calling

### Considerations before starting
Please, check the config file in the config folder and change the appropriate parts to use it with your own data, for example changing the name of the working directory or the name of the samples. The samples must be in the data/samples folder, and the genome, in the data folder.

See that along all the work we will be using this {} keys. In bash, it allows to exeute independently the thing inside that are separated by ','. For example if we have {a,b}.vcf it means a.vcf and b.vcf. 

Also the code below is written so it works with the samples in this repository, if you use other data you need to change the name in some parts of the cose. For example if you see A_1_20, and you use other samples, when you copy the code you need to change the names to use the name of your samples.

If a chromosme length is higher than 512Mb gatk will not work. You could divide the chromosomes to run properly gatk. You can check the length of the chromosmes using:
```
awk ' $2>=512000000 {print "One or more chromosmoes are bigger that 512MB" }' data/GCA_904849725.1_MorexV3_pseudomolecules.chrnames.fna.fai
```

The steps to do out variant calling are:

  1- Preparation of our data.

  2- Quality Control of our samples.

  3- Mapping against the reference genome.

  4- Marking duplicates and treating the readgroups.

  5- Perform the variant calling with freebayes and gatk.

  6- Perform the variant calling with deepvariantl.

  7- Filter vcf files


The whole code can be ran using this, but ius recommended to do it step by step:
```
snakemake -p -c4
```

  ### 1-Preparation of our data.

 
It is possible that the data is compressed (as in our case), so for that we make a rule to descompress our archives. We will use:
  ```
  snakemake data/samples/A_1_20_{1,2}.fastq -c4
  ```
Then we will get the fastq for each of our samples. Fastq is a text format to store biological sequences,obtained from the sequencing, and their quality scores. 


### 2- Quality Control of our samples.


To get the quality of each of the fastq you have to run this code:

```
snakemake data/samples/A_1_20_{1,2}_fastqc.{html,zip} -c4
```
Quality control is important because it give us specific information about the quality of the fastq to study, and can help us deciding if we have to do a preprocessing to our sample, like removing the adapters used in the sequencing.


### 3- Mapping against the reference genome.

The first step is to index the reference genome. To do that you can run this code: 

```
bwa index -a bwtsw data/genome.fasta
```
If you see the snakefile, you will see that is not necessary to execute the code above because executing the code below it will call the rule that indexes the reference genome.

This is a step that can take a lot of time (from hours to days depending of the genome and the computational capacity). You can estimate the time is going to take because the terminal will show you so often how many characters have it studied so far.

Then we can proceed with the mapping of our samples against the reference genome :
```
snakemake results/mapped/A_1_20.bam -c4
```
The BAM file is the compressed binary version of a SAM file, that is the format to save the mapping of the sequences.
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

Now we are going to mark the duplicated reads caused during the sequenciacing step, for example during the PCR. To do this we use:

```
snakemake results/markdup/A_1_20_nodup.{bam,txt} -c4
```

### 5- Perform the variant calling. 

Finally we can perform the variant calling, to do that we have to run :
```
snakemake calls/bayes/call_bayes.vcf -c4
```
To get the variants with freebayes :
```
snakemake calls/gatk/call_gatk.vcf -c4
```
To get the variants using GATK :

This is also a process that take some time.

### 6-Perform the variant calling with deepvariant.

To use deepvariant we have to run docker and pull an image from it. If docker is not installed use:
```
sudo apt -y update
sudo apt-get -y install docker.io
```
Then we can pull the image with:

```
BIN_VERSION="1.6.0"
sudo docker pull google/deepvariant:"${BIN_VERSION}"
```
Finally to get the deepvariant vcf you must do:

```
snakemake calls/deepvariant/call_deepvariant.vcf.gz -c4
```
See that the deepvariant algorithm give us the vcf file compressed.

### 7- Filter vcf files.

Now we are going to filter our vcf files to divide it in the ones with snps and the ones with indels. 

To do this we will use bcftools, but first there are some previous steps we must do.

First we need to compress our files into an bgzip format, that compress the file in blocks and allow indexes to be built against the compressed file.The deepvariant file is already compressed.

To do that we use:
```
snakemake calls/{bayes/call_bayes,gatk/call_gatk}_sorted.vcf.gz -c4

```

It could be that out vcf files are not well sorted, so we must sort them using:

```
snakemake calls/{bayes/call_bayes,gatk/call_gatk,deepvariant/call_deepvariant}_sorted.vcf.gz -c4

```

Then we need to index our vcf files to use bcftools with the following code:

```
snakemake calls/{bayes/call_bayes,gatk/call_gatk,deepvariant/call_deepvariant}_sorted.vcf.gz.csi -c4

```

And finally we use the following code to get the vcf files:
```
snakemake calls/{bayes/call_bayes_{snp,indel},gatk/call_gatk_{snp,indel},deepvariant/call_deepvariant_{snp,indel}}.vcf -c4
```
## Results 

Using VariantQC, provided by GATK, we can see the differences in the variant calling of the three algorithms used:
![dag-plot-included](https://github.com/iniverman/barley-variant-calling/blob/main/call_quality.png)

## Workflow distribution

We can also check the steps taken by the program to perform the calling: 


```
snakemake --dag  -n  | dot -Tsvg > dag.svg

```

![dag-plot-included](https://github.com/iniverman/barley-variant-calling/blob/main/dag.svg)







