
rule kb_ref:
    """
    Make reference for kallisto, including introns.
    """
    input:
        genome_fa=config['ref']['sequence'],
        gtf=config['ref']['annotation']
    output:
        index="results/kb_ref/index.idx",
        t2g="results/kb_ref/t2g.txt",
        f1="results/kb_ref/cdna.fa",
        f2="results/kb_ref/intron.fa",
        c1="results/kb_ref/cdna_t2c.txt",
        c2="results/kb_ref/intron_t2c.txt",
    benchmark:
        "benchmarks/kb_ref/bench.txt"
    params:
    threads: 32
    resources:
        mem_gb=96,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    envmodules:
        config['modules']['kb']
    shell:
        """
        kb ref --overwrite \
        -i {output.index} \
        -g {output.t2g} \
        -f1 {output.f1} -f2 {output.f2} \
        -c1 {output.c1} -c2 {output.c2} \
        --workflow lamanno \
        -t {threads} \
        {input.genome_fa} {input.gtf}
        """

rule kb_count:
    """
    Run kb count.
    """
    input:
        fqs=expand("results/renamed_fastqs/{{sample}}_R{read}.fastq.gz", read=["1","2"]),
        idx="results/kb_ref/index.idx",
        t2g="results/kb_ref/t2g.txt",
        c1="results/kb_ref/cdna_t2c.txt",
        c2="results/kb_ref/intron_t2c.txt"
    output:
        expand("results/kb_count/{{sample}}/{fn}", fn=['flens.txt', 'index.saved', 'inspect.json', 'kb_info.json', 'matrix.ec', 'output.bus', 'run_info.json', 'transcripts.txt']),
        expand("results/kb_count/{{sample}}/counts_unfiltered/{fn}", fn=['cells_x_tcc.ambiguous.mtx', 'cells_x_tcc.barcodes.txt', 'cells_x_tcc.ec.txt', 'cells_x_tcc.mature.mtx', 'cells_x_tcc.nascent.mtx', 'cells_x_tcc.total.mtx', 'cells_x_tcc.nucleus.mtx', 'cells_x_tcc.cell.mtx']),
    benchmark:
        "benchmarks/kb_count/{sample}.txt"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        genome_fa=config['ref']['sequence'],
        gtf=config['ref']['annotation']
    threads: 32
    resources:
        mem_gb=96,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    envmodules:
        config['modules']['kb']
    shell:
        """
        kb count --overwrite \
                -w NONE \
                -m {resources.mem_gb}G  \
                -t {threads} \
                -i {input.idx} \
                -g {input.t2g} \
                --workflow nac \
                -x "STORM-SEQ" \
                --tcc \
                --sum 'total' \
                -o {params.outdir} \
                -c1 {input.c1} -c2 {input.c2} \
                {input.fqs}
        """



rule kallisto_quant_tcc:
    """
    Run kallisto quant_tcc.
    """
    input:
        t2g="results/kb_ref/t2g.txt",
        flen="results/kb_count/{sample}/flens.txt",
        mtx="results/kb_count/{sample}/counts_unfiltered/cells_x_tcc.{type}.mtx",
        index="results/kb_ref/index.idx",
        ec="results/kb_count/{sample}/counts_unfiltered/cells_x_tcc.ec.txt",
    output:
         expand("results/kallisto_quant_tcc/{{sample}}/{{type}}/{fn}", fn=['genes.txt','matrix.abundance.gene.tpm.mtx','matrix.abundance.tpm.mtx','matrix.abundance.gene.mtx','matrix.abundance.mtx','transcripts.txt','abundance.gene_1.tsv','abundance_1.h5','abundance_1.tsv'])
    benchmark:
        "benchmarks/kallisto_quant_tcc/{sample}.{type}.txt"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0])
    threads: 1
    resources:
        mem_gb=96,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    envmodules:
        config['modules']['kallisto']
    shell:
        """
        kallisto quant-tcc \
                -i {input.index} \
                --genemap {input.t2g} \
                -f {input.flen} \
                -t {threads}  \
                -o {params.outdir} \
                --matrix-to-files \
                --ec-file {input.ec} \
                {input.mtx}

        """


