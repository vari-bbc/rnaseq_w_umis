rule make_erna_gtf:
    """
    Make eRNA GTF file.
    """
    input:
        genes = config["ref"]["annotation"],
        fantom5_enhancers = config["fantom5_enhancers"],
        genome = config["ref"]["fai"]
    output:
        genes_plus_2kb = "results/make_erna_gtf/genes_plus_2kb.gtf",
        nongenic_enh = "results/make_erna_gtf/fantom5_nongenic_enhancers.bed",
        nongenic_enh_genepred = "results/make_erna_gtf/fantom5_nongenic_enhancers.genepred",
        nongenic_enh_gtf = "results/make_erna_gtf/fantom5_nongenic_enhancers.gtf",
        genes_and_enhancers = "results/make_erna_gtf/genes_and_enhancers.gtf"
    params:
    threads: 16
    resources:
        mem_gb = 180,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    log:
    benchmark:
        "benchmarks/make_erna_gtf/bench.txt"
    envmodules:
        config["modules"]["bedtools"],
        config["modules"]["ucsc_tools"],
    shell:
       """
       bedtools slop -i {input.genes} -g {input.genome} -b 2000 > {output.genes_plus_2kb}
       bedtools intersect -a {input.fantom5_enhancers} -b {output.genes_plus_2kb} -v > {output.nongenic_enh}
       bedToGenePred {output.nongenic_enh} {output.nongenic_enh_genepred}
       genePredToGtf file {output.nongenic_enh_genepred} {output.nongenic_enh_gtf}
       cat {input.genes} {output.nongenic_enh_gtf} > {output.genes_and_enhancers}
       """

rule eRNAs_star_idx:
    input:
        genome_fa=config['ref']['sequence'],
        genes_gtf="results/make_erna_gtf/genes_and_enhancers.gtf"
    output:
        "results/eRNAs_star_idx/Log.out",
        "results/eRNAs_star_idx/SA",
        "results/eRNAs_star_idx/SAindex"
    log:
    benchmark:
        "benchmarks/eRNAs_star_idx/bench.txt"
    params:
        sjdb_overhang=100,
        outpref="results/eRNAs_star_idx/",
        ram=290000000000
    threads:16
    resources:
        mem_gb=300,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules:
        config["modules"]["star"]
    shell:
        """
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --outFileNamePrefix {params.outpref} \
            --genomeDir {params.outpref} \
            --genomeFastaFiles {input.genome_fa} \
            --sjdbGTFfile {input.genes_gtf} \
            --sjdbOverhang {params.sjdb_overhang} \
            --limitGenomeGenerateRAM {params.ram}

        """

rule STARsolo_eRNAs:
    """
    Run STARsolo.
    """
    input:
        unpack(get_star_solo_input),
        idx="results/eRNAs_star_idx/Log.out"
    output:
        star = multiext("results/STARsolo_eRNAs/{sample}.",
                 "Aligned.sortedByCoord.out.bam",
                 "Aligned.sortedByCoord.out.bam.bai",
                 "Log.final.out",
                 "Log.out",
                 "SJ.out.tab"),
        raw = expand("results/STARsolo_eRNAs/{{sample}}.Solo.out/{gene}/raw/{file}", gene=['Gene','GeneFull'], file=['matrix.mtx.gz','barcodes.tsv.gz','features.tsv.gz']),
    params:
        index = lambda wildcards, input: os.path.dirname(input.idx),
        outprefix = "results/STARsolo_eRNAs/{sample}.",
        tech_params = get_star_solo_params,
        # adjust R1 and R2 input order depending on tech.
        fqs = lambda wildcards, input: input.fq2 + ' ' + input.fq1 if(config["starsolo_tech"] in ["10x_v1", "10x_v2", "10x_v3", "cellseq192"]) else input.fq1 + ' ' + input.fq2
    threads: 16
    resources:
        mem_gb = 180,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    priority: 1
    log:
    benchmark:
        "benchmarks/STARSolo_eRNAs/{sample}.txt"
    envmodules:
        config["modules"]["star"],
        config["modules"]["samtools"],
        config["modules"]["pigz"]
    shell:
       """
       STAR  \
       --runThreadN {threads} \
       --limitBAMsortRAM 137438953472 \
       --genomeDir {params.index} \
       --readFilesIn {params.fqs} \
       --outSAMattributes NH HI AS nM CB UB \
       --readFilesCommand zcat  \
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix {params.outprefix} \
       {params.tech_params}

       samtools index {params.outprefix}Aligned.sortedByCoord.out.bam

       pigz -p {threads} {params.outprefix}Solo.out/Gene/raw/*
       pigz -p {threads} {params.outprefix}Solo.out/GeneFull/raw/*
       """

