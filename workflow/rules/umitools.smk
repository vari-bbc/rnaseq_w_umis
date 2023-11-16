
rule umitools_extract:
    """
    Extract UMIs.
    """
    input:
        fqs=expand("results/renamed_fastqs/{{sample}}_R{read}.fastq.gz", read=['1','2']),
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
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
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
    output:
        bam = "results/star/{stage}/{sample}.Aligned.sortedByCoord.out.bam",
        bai = "results/star/{stage}/{sample}.Aligned.sortedByCoord.out.bam.bai",
        counts = "results/star/{stage}/{sample}.ReadsPerGene.out.tab",
        flagstat = "results/star/{stage}/{sample}.Aligned.sortedByCoord.out.bam.flagstat",
    params:
        outprefix = "results/star/{stage}/{sample}.",
        index = config["ref"]["index"],
    benchmark:
        "benchmarks/star/{stage}/{sample}.txt"
    envmodules:
        config['modules']['star'],
        config['modules']['samtools']
    threads: 8
    resources:
        nodes =   1,
        mem_gb =  120,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {input.fqs} \
        --twopassMode Basic \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix {params.outprefix} \
        --quantMode GeneCounts \
        --outStd Log

        samtools index -@ {threads} {output.bam}
        samtools flagstat -@ {threads} {output.bam} > {output.flagstat}
        """


rule filt_raw_bams:
    """
    Filter undeduped BAMs.
    """
    input:
        bam = "results/star/raw/{sample}.Aligned.sortedByCoord.out.bam",
        bai = "results/star/raw/{sample}.Aligned.sortedByCoord.out.bam.bai",
    output:
        bam = "results/filt_raw_bams/{sample}.Aligned.sortedByCoord.out.bam",
        bai = "results/filt_raw_bams/{sample}.Aligned.sortedByCoord.out.bam.bai",
        flagstat = "results/filt_raw_bams/{sample}.Aligned.sortedByCoord.out.bam.flagstat",
    benchmark:
        "benchmarks/filt_raw_bams/{sample}.txt"
    params:
    threads: 8
    resources:
        mem_gb=72,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    envmodules:
        config['modules']['samtools']
    shell:
        """

        samtools view -@ {threads} -F 256 -f 2  -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        samtools flagstat -@ {threads} {output.bam} > {output.flagstat}

        """

rule umitools_dedup:
    """
    Extract UMIs.
    """
    input:
        bam = "results/filt_raw_bams/{sample}.Aligned.sortedByCoord.out.bam",
        bai = "results/filt_raw_bams/{sample}.Aligned.sortedByCoord.out.bam.bai",
    output:
        bam = "results/umitools_dedup/{sample}.Aligned.sortedByCoord.out.bam",
        bai = "results/umitools_dedup/{sample}.Aligned.sortedByCoord.out.bam.bai",
    benchmark:
        "benchmarks/umitools_dedup/{sample}.txt"
    params:
        sam="results/umitools_dedup/{sample}.Aligned.sortedByCoord.out.sam",
        tmp="/tmp",
        out_pref="results/umitools_dedup/{sample}"
    threads: 8
    resources:
        mem_gb=72,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
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
                --out-sam -S {params.sam}

        samtools view -@ {threads} -o {output.bam} {params.sam}

        samtools index -@ {threads} {output.bam}
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
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
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
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
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
