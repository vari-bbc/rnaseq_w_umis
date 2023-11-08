import pandas as pd
import os
from shutil import which
from snakemake.utils import validate, min_version
from pathlib import Path
##### set minimum snakemake version #####
min_version("7.25.0")

configfile: "bin/config.yaml"

samples = pd.read_table("bin/units.tsv")

rule all:
    input:
        expand("results/cogent_analyze/{sample}/", sample=pd.unique(samples['sample'])) if config['run_cogent_ap'] else [],
        #expand("results/htseq_count/{sample}.{strandedness}.out", sample=pd.unique(samples['sample']), strandedness=['yes','no','reverse']),
        "results/multiqc/multiqc_report.html",
        expand("results/salmon/{sample}/logs/salmon_quant.log", sample=pd.unique(samples['sample'])),
        "results/SummarizedExperiment/sce.rds"


rule multiqc:
    input:
        expand("results/fastqc/{sample}_R{read}_fastqc.html", sample=pd.unique(samples['sample']), read=['1','2']),
        expand("results/trim_galore/{sample}_R{read}_val_{read}_fastqc.html", sample=pd.unique(samples['sample']), read=["1","2"]),
        expand("results/fastq_screen/{sample}_R{read}_screen.html", sample=pd.unique(samples['sample']), read=['1','2']),
        expand("results/star/{stage}/{sample}.Aligned.sortedByCoord.out.bam.bai", stage=['raw','deduped'], sample=pd.unique(samples['sample'])),
        expand("results/salmon/{sample}/logs/salmon_quant.log", sample=pd.unique(samples['sample'])),
    output:
        "results/multiqc/multiqc_report.html",
    benchmark:
        "benchmarks/multiqc/multiqc.txt"
    params:
        workdir="results/multiqc/",
        dirs=lambda wildcards, input: ' '.join(pd.unique(['/'.join(Path(x).parts[0:2]) for x in input])),
        outfile="multiqc_report"
    envmodules:
        config['modules']['multiqc']
    threads: 4
    resources:
        mem_gb=100,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        multiqc \
        --force \
        --config bin/multiqc_config.yaml \
        --outdir {params.workdir} \
        --filename {params.outfile} \
        {params.dirs}

        """

def get_orig_fastq(wildcards):
    if wildcards.read == "R1":
            fastq = expand("raw_data/{fq}", fq = samples[samples["sample"] == wildcards.sample]["fq1"].values)
    elif wildcards.read == "R2":
            fastq = expand("raw_data/{fq}", fq = samples[samples["sample"] == wildcards.sample]["fq2"].values)
    return fastq 

rule rename_fastqs:
    """
    Rename fastqs by biologically meaningful name. Concatenate different runs of same library.
    """
    input:
        get_orig_fastq
    output:
        "results/renamed_data/{sample}_{read}.fastq.gz"
    benchmark:
        "benchmarks/rename_fastqs/{sample}_{read}.txt"
    params:
        num_input=lambda wildcards, input: len(input),
        input=lambda wildcards, input: ["'" + x + "'" for x in input]
    threads: 1
    resources:
        mem_gb=8,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules:
    shell:
        """
        if [ {params.num_input} -gt 1 ]
        then
            cat {params.input} > {output}
        else
            ln -sr {params.input} {output}
        fi

        """

rule fastqc:
    """
    Run fastqc on raw_data/ files.
    """
    input:
        "results/renamed_data/{fq_pref}.fastq.gz"
    output:
        html="results/fastqc/{fq_pref}_fastqc.html",
        zip="results/fastqc/{fq_pref}_fastqc.zip"
    params:
        outdir="results/fastqc/"
    benchmark:
        "benchmarks/fastqc/{fq_pref}.txt"
    envmodules:
        config['modules']['fastqc']
    threads: 1
    resources:
        mem_gb = 32,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        fastqc --outdir {params.outdir} {input}
        """

rule fastq_screen:
    """
    Run fastq_screen to detect any contamination from other species or excessive rRNA.
    """
    input:
        "results/renamed_data/{fq_pref}.fastq.gz"
    output:
        html = "results/fastq_screen/{fq_pref}_screen.html",
        txt = "results/fastq_screen/{fq_pref}_screen.txt",
    params:
    benchmark:
        "benchmarks/fastq_screen/{fq_pref}.txt"
    envmodules:
        config['modules']['fastq_screen']
    threads: 8
    resources:
        mem_gb = 32,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        fastq_screen --threads {threads} --outdir results/fastq_screen/ {input}
        """

rule umitools_extract:
    """
    Extract UMIs.
    """
    input:
        fqs=expand("results/renamed_data/{{sample}}_R{read}.fastq.gz", read=['1','2']),
    output:
        fqs=expand("results/umitools_extract/{{sample}}_R{read}.fastq.gz", read=['1','2']),
        log="results/umitools_extract/{sample}.log.txt"
    benchmark:
        "benchmarks/umitools_extract/{sample}.txt"
    params:
        uncompressed_fqs=lambda wildcards, output: [os.path.splitext(x)[0]  for x in output.fqs]
    threads: 8
    resources:
        mem_gb=72,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules:
        config['modules']['umi_tools'],
        config['modules']['pigz'],
    shell:
        """
        umi_tools extract \
                --bc-pattern=NNNNNNNN \
                --log={output.log} \
                -I {input.fqs[1]} \
                --read2-in={input.fqs[0]} \
                --stdout={params.uncompressed_fqs[1]} \
                --read2-out={params.uncompressed_fqs[0]}

        pigz -p {threads} {params.uncompressed_fqs[1]}
        pigz -p {threads} {params.uncompressed_fqs[0]}
        """

rule trim_galore_PE:
    input:
        expand("results/umitools_extract/{{sample}}_R{read}.fastq.gz", read=['1','2'])
    output:
        expand("results/trim_galore/{{sample}}_R{read}_val_{read}{suffix}", read=['1', '2'], suffix=['.fq.gz','_fastqc.html','_fastqc.zip']),
        expand("results/trim_galore/{{sample}}_R{read}.fastq.gz_trimming_report.txt", read=['1', '2']),
    params:
        outdir="results/trim_galore/",
        clip_R2=config['clip_r2']
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    envmodules:
        config['modules']['trim_galore']
    threads: 4
    resources:
        nodes =   1,
        mem_gb =  32,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        trim_galore \
        --paired \
        {input} \
        --output_dir {params.outdir} \
        --cores {threads} \
        {params.clip_R2} \
        -q 20 \
        --fastqc
        """


rule STAR:
    input:
        fqs=lambda wildcards: expand("results/trim_galore/{{sample}}_R{read}_val_{read}.fq.gz", read=['1','2']) if wildcards.stage == "raw" else expand("results/deduped_bam_to_fq/{{sample}}_dedup_R{read}.fq.gz", read=['1','2']),
        index = config["ref"]["index"],
    output:
        bam = "results/star/{stage}/{sample}.Aligned.sortedByCoord.out.bam",
        bai = "results/star/{stage}/{sample}.Aligned.sortedByCoord.out.bam.bai",
    params:
        outprefix = "results/star/{stage}/{sample}."
    benchmark:
        "benchmarks/star/{stage}/{sample}.txt"
    envmodules:
        config['modules']['star'],
        config['modules']['samtools']
    threads: 8
    resources:
        nodes =   1,
        mem_gb =  120,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --genomeDir {input.index} \
        --readFilesIn {input.fqs} \
        --twopassMode Basic \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix {params.outprefix} \
        --quantMode GeneCounts \
        --outStd Log 

        samtools index -@ {threads} {output.bam}
        """

rule umitools_dedup:
    """
    Extract UMIs.
    """
    input:
        bam = "results/star/raw/{sample}.Aligned.sortedByCoord.out.bam",
        bai = "results/star/raw/{sample}.Aligned.sortedByCoord.out.bam.bai",
    output:
        bam = "results/umitools_dedup/{sample}.Aligned.sortedByCoord.out.bam",
        bai = "results/umitools_dedup/{sample}.Aligned.sortedByCoord.out.bam.bai",
    benchmark:
        "benchmarks/umitools_dedup/{sample}.txt"
    params:
        tmp="/tmp",
        sam=lambda wildcards, output: os.path.splitext(output.bam)[0] + '.sam',
        out_pref="results/umitools_dedup/{sample}"
    threads: 8
    resources:
        mem_gb=72,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules:
        config['modules']['umi_tools'],
        config['modules']['samtools']
    shell:
        """
        umi_tools dedup \
                --paired \
                --temp-dir={params.tmp} \
                -I {input.bam} \
                --output-stats={params.out_pref} \
                --multimapping-detection-method="NH" \
                --out-sam \
                -S {params.sam}
        
        samtools view -@ {threads} -o {output.bam} {params.sam}
        samtools index -@ {threads} {output.bam}
        
        rm {params.sam}
        """

rule htseq_count:
    """
    htseq count.
    """
    input:
        bam = "results/umitools_dedup/{sample}.Aligned.sortedByCoord.out.bam",
        bai = "results/umitools_dedup/{sample}.Aligned.sortedByCoord.out.bam.bai",
        gtf=config['ref']['annotation']
    output:
        "results/htseq_count/{sample}.{strandedness}.out"
    benchmark:
        "benchmarks/htseq_count/{sample}.{strandedness}.txt"
    params:
        tmp="/tmp",
    threads: 1
    resources:
        mem_gb=72,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules:
        config['modules']['htseq'],
    shell:
        """
        htseq-count -f 'bam' -r 'pos' \
                --stranded={wildcards.strandedness} \
                --type='exon' \
                --idattr='gene_id' \
                --mode='union' \
                --nonunique='none' \
                {input.bam} \
                {input.gtf} > {output}
        """

rule deduped_bam_to_fq:
    """
    Convert deduped bams to fastq.
    """
    input:
        bam = "results/umitools_dedup/{sample}.Aligned.sortedByCoord.out.bam",
        bai = "results/umitools_dedup/{sample}.Aligned.sortedByCoord.out.bam.bai",
    output:
        fq1 = "results/deduped_bam_to_fq/{sample}_dedup_R1.fq.gz",
        fq2 = "results/deduped_bam_to_fq/{sample}_dedup_R2.fq.gz",
    benchmark:
        "benchmarks/deduped_bam_to_fq/{sample}.txt"
    params:
        uncomp_fq1 = lambda wildcards, output: os.path.splitext(output.fq1)[0],
        uncomp_fq2 = lambda wildcards, output: os.path.splitext(output.fq2)[0],
        tmp="/tmp",
    threads: 8
    resources:
        mem_gb=72,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules:
        config['modules']['bedtools'],
        config['modules']['pigz'],
    shell:
        """
        bedtools bamtofastq -i {input.bam} -fq {params.uncomp_fq1} -fq2 {params.uncomp_fq2}

        pigz -p {threads} {params.uncomp_fq1}
        pigz -p {threads} {params.uncomp_fq2}
        """


rule make_salmon_idx:
    input:
        transcripts = [config['salmon']['transcripts'], config['salmon']['spikeins']] if config['salmon']['spikeins'] else config['salmon']['transcripts'],
        decoys = config['salmon']['decoys'],
    output:
        index_files = expand("results/make_salmon_idx/{{ref_id}}/{fn}", fn=['ctable.bin','ctg_offsets.bin','pos.bin','mphf.bin','rank.bin','refseq.bin','seq.bin','ref_indexing.log']),
        decoys_txt = "results/make_salmon_idx/{ref_id}.decoys.txt",
        tx_fasta = "results/make_salmon_idx/{ref_id}.fa",
    params:
        out_idx="results/make_salmon_idx/{ref_id}"
    benchmark:
        "benchmarks/make_salmon_idx/{ref_id}.txt"
    envmodules:
        config['modules']['salmon']
    threads: 32
    resources:
        nodes =   1,
        mem_gb =  120,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        grep -Po '^>\S+' {input.decoys} | perl -npe 's/^>//' > {output.decoys_txt}
        cat {input.transcripts} {input.decoys} > {output.tx_fasta}

        salmon index -t {output.tx_fasta} -i {params.out_idx} --gencode -d {output.decoys_txt} -p {threads}
        """

rule salmon:
    input:
        fqs=expand("results/deduped_bam_to_fq/{{sample}}_dedup_R{read}.fq.gz", read=["1","2"]),
        index=expand("results/make_salmon_idx/{ref_id}/{fn}", ref_id=config['salmon']['index_id'], fn=['ctable.bin','ctg_offsets.bin','pos.bin','mphf.bin','rank.bin','refseq.bin','seq.bin','ref_indexing.log'])
    output:
        expand("results/salmon/{{sample}}/{file}", file=["libParams/flenDist.txt","aux_info/meta_info.json","quant.sf","lib_format_counts.json","cmd_info.json","logs/salmon_quant.log"])
    params:
        outdir=directory("results/salmon/{sample}"),
        index_dir = lambda wildcards, input: '/'.join(Path(input.index[0]).parts[0:3])
    benchmark:
        "benchmarks/salmon/{sample}.txt"
    envmodules:
        config['modules']['salmon']
    threads: 8
    resources:
        nodes =   1,
        mem_gb =  120,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        salmon quant \
                -p {threads} \
                -l A \
                -i {params.index_dir} \
                -1 {input.fqs[0]} -2 {input.fqs[1]} \
                --validateMappings \
                -o {params.outdir}
        """


rule make_barcodes_file:
    input:
        "results/renamed_data/{sample}_R1.fastq.gz"
    output:
        "results/make_barcodes_file/{sample}.barcode"
    benchmark:
        "benchmarks/make_barcodes_file/{sample}.txt"
    params:
        sample=lambda wildcards: wildcards.sample
    threads: 1
    resources:
        mem_gb=48,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules:
    shell:
        """
        printf "Sample\tBarcode\n" > {output} 
        paste <(echo {params.sample}) <(zcat {input} | head -n 400000 | perl -lne '/(1:N:0:)(\S+$)/; print $2 if (($. % 4) == 1)' | sort | uniq -c | perl -npe 's/^\s+//; s/ /\\t/' | sort -k1,1nr | head -n1 | grep -Po '\S+$') >> {output}
        """

rule cogent_demux:
    input:
        fqs=expand("results/renamed_data/{{sample}}_R{read}.fastq.gz", read=['1','2']),
        barcode="results/make_barcodes_file/{sample}.barcode"
    output:
        directory("results/cogent_demux/{sample}/")
    benchmark:
        "benchmarks/cogent_demux/{sample}.txt"
    params:
        outdir="results/cogent_demux/{sample}",
        cogent_tech=config['cogent_ap_tech']
    threads: 8
    resources:
        mem_gb=72,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules:
        config['modules']['cogent_ap']
    shell:
        """
        # Remove the directory that snakemake autocreates or cogent errors out
        #rm -r {params.outdir}

        cogent demux \
        -i {input.fqs[0]} \
        -p {input.fqs[1]} \
        -t {params.cogent_tech} \
        -b {input.barcode} \
        -o {params.outdir} \
        -n {threads}
        """

rule cogent_analyze:
    input:
        "results/cogent_demux/{sample}/"
    output:
        directory("results/cogent_analyze/{sample}")
    benchmark:
        "benchmarks/cogent_analyze/{sample}.txt"
    params:
        outdir="results/cogent_analyze/{sample}",
        cogent_ref=config['cogent_ap_ref'],
        cogent_tech=config['cogent_ap_tech']
    threads: 8
    resources:
        mem_gb=72,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules:
        config['modules']['cogent_ap']
    shell:
        """
        # Remove the directory that snakemake autocreates or cogent errors out
        # rm -r {params.outdir}

        cogent analyze \
        -i {input} \
        -g {params.cogent_ref} \
        -t {params.cogent_tech} \
        -o {params.outdir} \
        --threads {threads}
        """

rule SummarizedExperiment:
    input:
        star = expand("results/star/deduped/{samples.sample}.ReadsPerGene.out.tab", samples=samples.itertuples()),
        salmon = expand("results/salmon/{samples.sample}/{file}", samples=samples.itertuples(), file=['lib_format_counts.json', 'quant.sf'])
    output:
        se="results/SummarizedExperiment/SummarizedExperiment.rds",
        sce="results/SummarizedExperiment/sce.rds",
        sizeFactors="results/SummarizedExperiment/DESeq2_sizeFactors_reciprocal.tsv"
    benchmark:
        "benchmarks/SummarizedExperiment/SummarizedExperiment.txt"
    params:
        gtf=config['ref']['annotation'],
        orgdb=config['orgdb']
    threads: 1
    envmodules:
        config['modules']['R']
    resources:
        mem_gb=64,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        Rscript --vanilla workflow/scripts/make_sce.R {params.gtf} {params.orgdb} {output.se} {output.sce} {output.sizeFactors}
        """
