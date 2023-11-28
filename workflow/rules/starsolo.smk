

rule add_fake_BC:
    """
    Add fake cell barcode for STARsolo.
    """
    input:
        fqs=expand("results/renamed_fastqs/{{sample}}_R{read}.fastq.gz", read=['1','2']),
    output:
        fqs=temp(expand("results/add_fake_BC/{{sample}}_R{read}.fastq.gz", read=['1','2'])),
    benchmark:
        "benchmarks/add_fake_BC/{sample}.txt"
    params:
        uncompressed_fqs=lambda wildcards, output: [os.path.splitext(x)[0]  for x in output.fqs]
    threads: 8
    resources:
        mem_gb=72,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    envmodules:
        config['modules']['pigz'],
    shell:
        """
        # Keep R1 the same
        ln -sr {input.fqs[0]} {output.fqs[0]}

        # Add A to the 5' end of R2; use for both the base and the BQ
        zcat {input.fqs[1]} | perl -npe 'if (($. % 4 == 2) || ($. % 4 == 0)){{$_ = "A".$_}}' > {params.uncompressed_fqs[1]}

        pigz -p {threads} {params.uncompressed_fqs[1]}
        """


def get_star_solo_params(wildcards):
    # modify parameters depending on kit/tech

    # In STORM, R2 begins with 8 bp UMI then 6 bp of linker/adapter
    # STARsolo requires a CB, so the input for STAR have 1 bp prepended to R2
    if(config["starsolo_tech"] == "STORM"):
        star_solo_params = """--soloStrand Reverse \
        --soloCellFilter None \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist None \
        --soloCBmatchWLtype Exact \
        --soloBarcodeMate 2 \
        --clip5pNbases 0 15 \
        --soloCBstart 1 \
        --soloCBlen 1 \
        --soloUMIstart 2 \
        --soloUMIlen 8 \
        --soloBarcodeReadLength 0 \
        --soloUMIdedup Exact \
        --soloFeatures Gene GeneFull"""

    return(star_solo_params)

def get_star_solo_input(wildcards):
    reads_dir = ""

    # Set the fastq file paths
    if(config["starsolo_tech"] in ["STORM"]):
        reads_dir = "results/add_fake_BC"
    else:
        reads_dir = "results/renamed_fastqs"

    star_solo_input = {}
    star_solo_input['fq1'] = reads_dir + "/{sample}_R1.fastq.gz".format(sample=wildcards.sample)
    star_solo_input['fq2'] = reads_dir + "/{sample}_R2.fastq.gz".format(sample=wildcards.sample)

    return(star_solo_input)

rule STARsolo:
    """
    Run STARsolo.
    """
    input:
        unpack(get_star_solo_input),
    output:
        star = multiext("results/STARsolo/{sample}.",
                 "Aligned.sortedByCoord.out.bam",
                 "Aligned.sortedByCoord.out.bam.bai",
                 "Log.final.out",
                 "Log.out",
                 "SJ.out.tab"),
        raw = expand("results/STARsolo/{{sample}}.Solo.out/{gene}/raw/{file}", gene=['Gene','GeneFull'], file=['matrix.mtx.gz','barcodes.tsv.gz','features.tsv.gz']),
    params:
        index = config["ref"]["index"],
        outprefix = "results/STARsolo/{sample}.",
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
        "benchmarks/STARSolo/{sample}.txt"
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


