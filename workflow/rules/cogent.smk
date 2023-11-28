
rule make_barcodes_file:
    input:
        "results/renamed_fastqs/{sample}_R1.fastq.gz"
    output:
        "results/make_barcodes_file/{sample}.barcode"
    benchmark:
        "benchmarks/make_barcodes_file/{sample}.txt"
    params:
        sample=lambda wildcards: wildcards.sample
    threads: 1
    resources:
        mem_gb=48,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    envmodules:
    shell:
        """
        printf "Sample\tBarcode\n" > {output}
        paste <(echo {params.sample}) <(zcat {input} | head -n 400000 | perl -lne '/(1:N:0:)(\S+$)/; print $2 if (($. % 4) == 1)' | sort | uniq -c | perl -npe 's/^\s+//; s/ /\\t/' | sort -k1,1nr | head -n1 | grep -Po '\S+$') >> {output}
        """

rule cogent_demux:
    input:
        fqs=expand("results/renamed_fastqs/{{sample}}_R{read}.fastq.gz", read=['1','2']),
        barcode="results/make_barcodes_file/{sample}.barcode"
    output:
        temp(directory("results/cogent_demux/{sample}/"))
    benchmark:
        "benchmarks/cogent_demux/{sample}.txt"
    params:
        outdir="results/cogent_demux/{sample}",
        cogent_tech=config['cogent_ap_tech']
    threads: 8
    resources:
        mem_gb=72,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
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
        expand("results/cogent_analyze/{{sample}}/{{sample}}{file_suff}", file_suff=["_umi_uss_genematrix.csv","_analyzer.log","_stats.csv"]),
        "results/cogent_analyze/{sample}/gene_info.csv",
        temp(directory("results/cogent_analyze/{sample}/bam")),
        directory("results/cogent_analyze/{sample}/extras"),
        directory("results/cogent_analyze/{sample}/work"),
    benchmark:
        "benchmarks/cogent_analyze/{sample}.txt"
    params:
        outdir="results/cogent_analyze/{sample}",
        cogent_ref=config['cogent_ap_ref'],
        cogent_tech=config['cogent_ap_tech']
    threads: 8
    resources:
        mem_gb=72,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    envmodules:
        config['modules']['cogent_ap']
    shell:
        """
        # Remove the directory that snakemake autocreates or cogent errors out
        rm -r {params.outdir}

        cogent analyze \
        -i {input} \
        -g {params.cogent_ref} \
        -t {params.cogent_tech} \
        -o {params.outdir} \
        --threads {threads}
        """
