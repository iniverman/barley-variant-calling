configfile: "config/config.yaml"

# This rule allow to run all the code bellow by giving as inputs the outputs 
# of some o the other rules

rule all:
    input: 
        expand("data/samples/{sample}_1_fastqc.html",sample=config["samples_name"]),
        expand("data/samples/{sample}_2_fastqc.html",sample=config["samples_name"]),
        expand("data/{genome_name}.bwt",genome_name=config["genome_name"]),
        expand("data/{genome_name}.fai",genome_name=config["genome_name"]),
        expand("calls/gatk/{vcf_gatk_name}.vcf",vcf_gatk_name=config["vcf_gatk_name"]),
        expand("calls/bayes/{vcf_bayes_name}.vcf",vcf_bayes_name=config["vcf_bayes_name"]),
        expand("calls/deepvariant/{vcf_deepvariant_name}.vcf.gz",vcf_deepvariant_name=config["vcf_deepvariant_name"]),
        expand("calls/bayes/{vcf_bayes_name}_snp.vcf",vcf_bayes_name=config["vcf_bayes_name"]),
        expand("calls/bayes/{vcf_bayes_name}_indel.vcf",vcf_bayes_name=config["vcf_bayes_name"]),
        expand("calls/gatk/{vcf_gatk_name}_snp.vcf",vcf_gatk_name=config["vcf_gatk_name"]),
        expand("calls/gatk/{vcf_gatk_name}_indel.vcf",vcf_gatk_name=config["vcf_gatk_name"])
        
        

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

#The next rule is to index the genome for bwa, as input it will take the
#genome and as output the indexes.
rule index_bwa_genome:
    input: 
        expand("data/{genome_name}",genome_name=config["genome_name"])
    output:
        expand("data/{genome_name}.bwt",genome_name=config["genome_name"])
    shell:
        "bwa index -a bwtsw {input}"
#The next rule is to do the same but with samtools, as input is the genome and 
#as output the indexes in .fai format.
rule index_samtools_genome:
    input: 
        expand("data/{genome_name}",genome_name=config["genome_name"])
    output:
        expand("data/{genome_name}.fai",genome_name=config["genome_name"])
    shell:
        "samtools faidx {input}"
# Map the samples against the reference genome. As input has the reference 
# genome and the couple of paired end samples and as output the BAM with the mapping.

rule bwa_map:
    input:
        index=expand("data/{genome_name}.bwt",genome_name=config["genome_name"]),
        ref=expand("data/{genome_name}",genome_name=config["genome_name"]),
        s1=expand("data/samples/{sample}_1.fastq",sample=config["samples_name"]),
        s2=expand("data/samples/{sample}_2.fastq",sample=config["samples_name"])
    output:
       temp("results/mapped/{sample}.bam")
    params:
        rg=r"-R '@RG\tID:Seq01a\tSM:sample_A_1_21\tPL:ILLUMINA'"
    threads: 16
    shell:
        "bwa mem {params.rg} -t {threads}  {input.ref} {input.s1} {input.s2} | samtools view -Sb - > {output}"

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

# Print the mapping depth of each chromosome in a subplot. As input it takes the previous csv
# and the output is the plot of subplots. 

rule plot_depth:
    input:
        expand("results/depth/{sample}_depth.csv",sample=config["samples_name"])
    output:
        "results/depth/plots/{sample}.svg"
    script:
        "scripts/plot-depth.py"     

# Mark duplicate reads of the mapping. The input is the BAM file sorted and as output
# is the BAM file marked and a txt with the metrics of the process.

rule mark_duplicates:
    input: 
        expand("results/sorted/{sample}_sorted.bam",sample=config["samples_name"])
    output:
        marked_bam="results/markdup/{sample}_nodup.bam",
        metrics="results/markdup/{sample}_nodup.txt"
    shell:
            "java -jar ../picard.jar MarkDuplicates -I {input} -O {output.marked_bam} -M {output.metrics}" 


#Index the marked BAM file, as input takes the BAM file and the output is the BAM file indexed.

rule samtools_index:
    input:
        marked_bam="results/markdup/{sample}_nodup.bam"
    output:
        "results/markdup/{sample}_nodup.bam.csi"
    shell:
        "samtools index -c {input.marked_bam} {output}"

# It does the variant calling using gatk. As input it takes the BAM file with the new 
# read groups and the reference genome and the output is the vcf file with the variant
# calling.        

rule haplotype_caller:
    input:
        index=expand("results/markdup/{sample}_nodup.bam.csi",sample=config["samples_name"]),
        bam=expand("results/markdup/{sample}_nodup.bam",sample=config["samples_name"]),
        ref=expand("data/{genome_name}",genome_name=config["genome_name"])
    output:
        expand("calls/gatk/{vcf_gatk_name}.vcf",vcf_gatk_name=config["vcf_gatk_name"])
    threads: 8 
    shell:
        " ../gatk-4.5.0.0/gatk --java-options  '-Xmx4g' " 
        "HaplotypeCaller -R {input.ref} -I {input.bam} -O {output} -ERC GVCF"

# It does the variant calling using freebayes. As input it takes the BAM file with the new 
# read groups and the reference genome and the output is the vcf file with the variant
# calling.

rule free_bayes:
    input:
        bam=expand("results/markdup/{sample}_nodup.bam",sample=config["samples_name"]),
        ref=expand("data/{genome_name}",genome_name=config["genome_name"])
    output:
        expand("calls/bayes/{vcf_bayes_name}.vcf",vcf_bayes_name=config["vcf_bayes_name"])
    shell:
        "freebayes -f {input.ref} {input.bam} --gvcf -g 1000 >{output}"

# This rule make the variant calling using deepvariant, it take as inputs the bam and 
#the bam indexed and the reference genome and the genome indexed. The output is the vcf file.
rule deepvariant:
    input:
        ref_index=expand("data/{genome_name}.fai",genome_name=config["genome_name"]),
        ref=expand("data/{genome_name}",genome_name=config["genome_name"]),
        bam=expand("results/markdup/{sample}_nodup.bam",sample=config["samples_name"]),
        bam_index=expand("results/markdup/{sample}_nodup.bam.csi",sample=config["samples_name"])
    output:
        expand("calls/deepvariant/{vcf_deepvariant_name}.vcf.gz",vcf_deepvariant_name=config["vcf_deepvariant_name"])
    shell:
        " mv  {input.ref_index} /data_deepvariant/"  
        "| mv -f {input.ref} /data_deepvariant |mv -f {input.bam_index} "
        "/data_deepvariant | mv -f {input.bam} /data_deepvariant |sudo docker run  "
        "-v /data_deepvariant:/input  "
        "-v /calls/deepvariant:/output "
        "google/deepvariant:1.6.0 "
        "/opt/deepvariant/bin/run_deepvariant  "
        "--model_type=WGS  "
        "--ref={input.ref}"
        "--reads={input.bam}"
        "--output_vcf=/output/output.vcf.gz  "
        "--output_gvcf=/output/output.g.vcf.gz  "
        "--intermediate_results_dir /output/intermediate_results_dir  "
        "--num_shards=1"
# Now we will compress the files in a bgzip file to use it with bcftools view. As input
# takes the vcf files and as outputo the compressed files.

rule compress:
    input: 
        bayes=expand("calls/bayes/{vcf_bayes_name}.vcf",vcf_bayes_name=config["vcf_bayes_name"]),
        gatk=expand("calls/gatk/{vcf_gatk_name}.vcf",vcf_gatk_name=config["vcf_gatk_name"])
    output: 
        gz1=expand("calls/bayes/{vcf_bayes_name}.vcf.gz",vcf_bayes_name=config["vcf_bayes_name"]),
        gz2=expand("calls/gatk/{vcf_gatk_name}.vcf.gz",vcf_gatk_name=config["vcf_gatk_name"])
    shell:
        "bgzip -c {input.bayes} > {output.gz1} | bgzip -c {input.gatk} > {output.gz2}"


# Next rule will sort the bayes file, to do then the indexes. As input it takes the bgzip compressed file and as output the 
# sorted vcf file

rule sort_bayes: 
    input:
        bayes=expand("calls/bayes/{vcf_bayes_name}.vcf.gz",vcf_bayes_name=config["vcf_bayes_name"]),
        
    output:
        s1=expand("calls/bayes/{vcf_bayes_name}_sorted.vcf.gz",vcf_bayes_name=config["vcf_bayes_name"]),
        
    shell: 
        "bcftools sort {input.bayes} -o {output.s1}"

# Now we can index the vf file with the following rule. As input it will take the sorted vcf
# and as output will creat a file with the indexes.

rule index:
    input:
        bayes=expand("calls/bayes/{vcf_bayes_name}_sorted.vcf.gz",vcf_bayes_name=config["vcf_bayes_name"]),
        gatk=expand("calls/gatk/{vcf_gatk_name}.vcf.gz",vcf_gatk_name=config["vcf_gatk_name"]),
        deepvariant=expand("calls/deepvariant/{vcf_deepvariant_name}.vcf.gz",vcf_deepvariant_name=config["vcf_deepvariant_name"])
    output:
        i1=expand("calls/bayes/{vcf_bayes_name}_sorted.vcf.gz.csi",vcf_bayes_name=config["vcf_bayes_name"]),
        i2=expand("calls/gatk/{vcf_gatk_name}.vcf.gz.csi",vcf_gatk_name=config["vcf_gatk_name"]),
        i3=expand("calls/deepvariant/{vcf_deepvariant_name}.vcf.gz.csi",vcf_deepvariant_name=config["vcf_deepvariant_name"])
    shell:
        "bcftools index {input.bayes} | bcftools index {input.gatk} | bcftools index {input.deepvariant}" 

# The next rule will divide the vcf in snp and indels. As output will take the compressed vcf
# and output the vcf divided.

rule divide_vcf: 
    input:
        i1=expand("calls/bayes/{vcf_bayes_name}_sorted.vcf.gz.csi",vcf_bayes_name=config["vcf_bayes_name"]),
        i2=expand("calls/gatk/{vcf_gatk_name}.vcf.gz.csi",vcf_gatk_name=config["vcf_gatk_name"]),
        i3=expand("calls/deepvariant/{vcf_deepvariant_name}.vcf.gz.csi",vcf_deepvariant_name=config["vcf_deepvariant_name"]),
        bayes=expand("calls/bayes/{vcf_bayes_name}_sorted.vcf.gz",vcf_bayes_name=config["vcf_bayes_name"]),
        gatk=expand("calls/gatk/{vcf_gatk_name}.vcf.gz",vcf_gatk_name=config["vcf_gatk_name"]),
        deepvariant=expand("calls/deepvariant/{vcf_deepvariant_name}.vcf.gz",vcf_deepvariant_name=config["vcf_deepvariant_name"])
    output:
        bayes_snp = expand("calls/bayes/{vcf_bayes_name}_snp.vcf",vcf_bayes_name=config["vcf_bayes_name"]),
        bayes_indel = expand("calls/bayes/{vcf_bayes_name}_indel.vcf",vcf_bayes_name=config["vcf_bayes_name"]),
        gatk_snp = expand("calls/gatk/{vcf_gatk_name}_snp.vcf",vcf_gatk_name=config["vcf_gatk_name"]),
        gatk_indel = expand("calls/gatk/{vcf_gatk_name}_indel.vcf",vcf_gatk_name=config["vcf_gatk_name"]),
        deepvariant_snp = expand("calls/deepvariant/{vcf_deepvariant_name}_snp.vcf",vcf_deepvariant_name=config["vcf_deepvariant_name"]),
        deepvariant_indel =expand("calls/deepvariant/{vcf_deepvariant_name}_indel.vcf",vcf_deepvariant_name=config["vcf_deepvariant_name"])
    shell: 
        "bcftools view -v snps {input.bayes} -Oz -o {output.bayes_snp} "
        "| bcftools view -v indels {input.bayes} -Oz -o {output.bayes_indel}"
        "| bcftools view -v snps {input.gatk} -Oz -o {output.gatk_snp}"
        "| bcftools view -v indels {input.gatk} -Oz -o {output.gatk_indel}"
        "| bcftools view -v snps {input.deepvariant} -Oz -o {output.deepvariant_snp}"
        "| bcftools view -v indels {input.deepvariant} -Oz -o {output.deepvariant_indel}"


#

