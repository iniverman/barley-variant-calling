configfile: "config/config.yaml"


rule all:
    input: 
        "data/samples/A_1_20_1_fastqc.html",
        "data/samples/A_1_20_2_fastqc.html",
        "calls/call_bayes.vcf"
        
        


rule unzip: 
    input: 
        sample1=expand("data/samples/{name}_1.fastq.gz",name=config["samples_name"]),
        sample2=expand("data/samples/{name}_2.fastq.gz",name=config["samples_name"])
    output:
        fq1=expand("data/samples/{name}_1.fastq",name=config["samples_name"]),
        fq2=expand("data/samples/{name}_2.fastq",name=config["samples_name"])
    shell:
        "gunzip -c {input.sample1} > {output.fq1} | gunzip -c {input.sample2} > {output.fq2}"


rule quality_control:
    input: 
        expand("data/samples/{name}_{unit}.fastq",name=config["samples_name"],unit=config["samples_unit"])
    output:
        html="data/samples/{name}_{unit}_fastqc.html",
        zip="data/samples/{name}_{unit}_fastqc.zip"
    shell:
        "fastqc {input}"

rule index_genome:
    input:
        genome= expand("/data/{genome}",genome=config["genome_name"]) 
    shell:
        "bwa index -a bwtsw {input.genome}"
        
rule bwa_map:
    input:
        genome=config["genome_name"],
        read1=expand("data/samples/{sample}_1.fastq",sample=config["samples_name"]),
        read2=expand("data/samples/{sample}_2.fastq",sample=config["samples_name"])
    output:
       temp("results/mapped/{sample}.bam")
    params:
        rg=r"-R '@RG\tID:Seq01a\tSM:sample_A_1_21\tPL:ILLUMINA'"
    threads: 16
    shell:
        "bwa mem {input.genome} {input.read1} {input.read2} | samtools view -Sb - > {output}"

rule sort_bam:
    input:
        expand("results/mapped/{sample}.bam",sample=config["samples_name"])
    output:
        "results/sorted/{sample}_sorted.bam"
    conda:
        "env/mapping.yaml"
    shell:
        "samtools sort -o {output} {input}"

rule depth_calc:
    input:
        expand("results/sorted/{sample}_sorted.bam",sample=config["samples_name"])
    output:
        "results/depth/{sample}_depth.csv"
    shell:
        "samtools depth  {input} > {output} "

rule mean:
    input: 
        expand("results/depth/{sample}_depth.csv",sample=config["samples_name"])
    output:
        "results/depth/{sample}_mean.csv"
    script:
        "scripts/new_column.py"

rule plot_depth:
    input:
        expand("results/depth/{sample}_depth.csv",sample=config["samples_name"])
    output:
        "results/depth/plots/{sample}.svg"
    script:
        "scripts/plot-depth.py"      

rule mark_duplicates:
    input: 
        expand("results/sorted/{sample}_sorted.bam",sample=config["samples_name"])
    output:
        marked_bam="results/markdup/{sample}_nodup.bam",
        metrics="results/markdup/{sample}_nodup.txt"
    shell:
        "java -jar /home/iveron/Descargas/picard.jar "
        "MarkDuplicates -I {input} -O {output.marked_bam} -M {output.metrics}"  

rule add_or_replace_RG:
    input:
        expand("results/markdup/{sample}_nodup.bam",sample=config["samples_name"])
    output:
        "results/markdup/{sample}_RG.bam"
    shell:
        "java -jar /home/iveron/Descargas/picard.jar AddOrReplaceReadGroups"
        " -I {input}  -O {output} -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 20 "

rule free_bayes:
    input:
        bam=expand("results/markdup/{sample}_RG.bam",sample=config["samples_name"]),
        ref=config["genome_name"]
    output:
        "calls/call_bayes.vcf"
    shell:
        "freebayes -f {input.ref} {input.bam} --gvcf -g 1000  > {output}"
