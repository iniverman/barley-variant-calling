configfile: "config/config.yaml"

# This rule allow to run all the code bellow by giving as inputs the outputs 
# of some o the other rules

rule all:
    input: 
        expand("data/samples/{sample}_1_fastqc.html",sample=config["samples_name"]),
        expand("data/samples/{sample}_2_fastqc.html",sample=config["samples_name"]),
        "calls/call_bayes.vcf",
        "calls/call_gatk.vcf"
        
        

# Decompress the files. Given as input the files compressed uses gunzip to 
# decompress them.

rule unzip: 
    input: 
        sample1=expand("data/samples/{name}_1.fastq.gz",name=config["samples_name"]),
        sample2=expand("data/samples/{name}_2.fastq.gz",name=config["samples_name"])
    output:
        fq1=expand("data/samples/{name}_1.fastq",name=config["samples_name"]),
        fq2=expand("data/samples/{name}_2.fastq",name=config["samples_name"])
    shell:
        "gunzip -c {input.sample1} > {output.fq1} | gunzip -c {input.sample2} > {output.fq2}"

# Create the fastqc for the samples. As input has the samples in fastq and output
# the html and zip with the quality control of the samples.

rule quality_control:
    input: 
        expand("data/samples/{name}_{unit}.fastq",name=config["samples_name"],unit=config["samples_unit"])
    output:
        html="data/samples/{name}_{unit}_fastqc.html",
        zip="data/samples/{name}_{unit}_fastqc.zip"
    shell:
        "fastqc {input}"

# Map the samples against the reference genome. As input has the reference 
# genome and the couple of paired end samples and as output the BAM with the mapping.

rule bwa_map:
    input:
        genome=config["genome_name"],
        s1=expand("data/samples/{sample}_1.fastq",sample=config["samples_name"]),
        s2=expand("data/samples/{sample}_2.fastq",sample=config["samples_name"])
    output:
       temp("results/mapped/{sample}.bam")
    params:
        rg=r"-R '@RG\tID:Seq01a\tSM:sample_A_1_21\tPL:ILLUMINA'"
    threads: 16
    shell:
        "bwa mem {input.genome} {input.s1} {input.s2} | samtools view -Sb - > {output}"

# Sort the BAM file get from the previous rule, as input it take the BAM file and as output 
# the BAM file sorted.      

rule sort_bam:
    input:
        expand("results/mapped/{sample}.bam",sample=config["samples_name"])
    output:
        "results/sorted/{sample}_sorted.bam"
    conda:
        "env/mapping.yaml"
    shell:
        "samtools sort -o {output} {input}"

# Calculate the depth of the mapping for each position in each chromosome. It has as input
# the sorted BAM file and as output the csv with the calculations.

rule depth_calc:
    input:
        expand("results/sorted/{sample}_sorted.bam",sample=config["samples_name"])
    output:
        "results/depth/{sample}_depth.csv"
    shell:
        "samtools depth  {input} > {output} "

#Print the mapping depth of each chromosome in a subplot. As input it takes the previous csv
# and the output is the plot of subplots. 

rule plot_depth:
    input:
        expand("results/depth/{sample}_depth.csv",sample=config["samples_name"])
    output:
        "results/depth/plots/{sample}.svg"
    script:
        "scripts/plot-depth.py"     

#Mark duplicate reads of the mapping. The input is the BAM file sorted and as output
# is the BAM file marked and a txt with the metrics of the process.

rule mark_duplicates:
    input: 
        expand("results/sorted/{sample}_sorted.bam",sample=config["samples_name"])
    output:
        marked_bam="results/markdup/{sample}_nodup.bam",
        metrics="results/markdup/{sample}_nodup.txt"
    shell:
        "java -jar /home/iveron/Descargas/picard.jar "
        "MarkDuplicates -I {input} -O {output.marked_bam} -M {output.metrics}" 

# Assign the reads to a new readgroup. As input takes the BAM fiole with the marked 
# duplicated and as output the new readgroup.

rule add_or_replace_RG:
    input:
        expand("results/markdup/{sample}_nodup.bam",sample=config["samples_name"])
    output:
        "results/markdup/{sample}_RG.bam"
    shell:
        "java -jar /home/iveron/Descargas/picard.jar AddOrReplaceReadGroups"
        " -I {input}  -O {output} -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 20 "

# It does the variant calling using fgatk. As input it takes the BAM file with the new 
# read groups and the reference genome and the output is the vcf file with the variant
#calling.        

rule haplotype_caller:
    input:
        bam=expand("results/markdup/{sample}_RG.bam",sample=config["samples_name"]),
        ref=config["genome_name"]
    output:
        "calls/call_gatk.vcf"
    threads: 8 
    shell:
        "/home/iveron/Descargas/gatk-4.5.0.0/gatk --java-options  '-Xmx4g' " 
        "HaplotypeCaller -R {input.ref} -I {input.bam} -O {output} -ERC GVCF" 

# It does the variant calling using freebayes. As input it takes the BAM file with the new 
# read groups and the reference genome and the output is the vcf file with the variant
#calling.

rule free_bayes:
    input:
        bam=expand("results/markdup/{sample}_RG.bam",sample=config["samples_name"]),
        ref=config["genome_name"]
    output:
        "calls/call_bayes.vcf"
    shell:
        "freebayes -f {input.ref} {input.bam} --gvcf -g 1000  > {output}"
