rule pull_genes:
    """
    Pull genes of interest from IsoSeq results
    """
    input:
        fasta=join(RESULTSDIR,"collapse","{sample}.hq.collapsed.fasta"), # Full fasta headers per isoform
        classification_filtered=join(RESULTSDIR,"classify","{sample}_classification.filtered_lite_classification.txt"),
        gff=join(RESULTSDIR,"classify","{sample}.sorted.filtered_lite.gff")
    output:
        classification=join(RESULTSDIR,"genes_of_interest","{sample}.{gene}.classification.filtered_lite_classification.txt"),
        gff=join(RESULTSDIR,"genes_of_interest","{sample}.{gene}.sorted.filtered_lite.gff"),
        tmp_list=temp(join(RESULTSDIR,"genes_of_interest","{sample}.{gene}.tmp")),
        id_list=temp(join(RESULTSDIR,"genes_of_interest","{sample}.{gene}.list")),
        fasta=join(RESULTSDIR,"genes_of_interest","{sample}.{gene}.fasta"),
        bam=temp(join(RESULTSDIR,"genes_of_interest","{sample}.{gene}.bam")),
        sorted_bam=join(RESULTSDIR,"genes_of_interest","{sample}.{gene}.sorted.bam")
    envmodules:
        TOOLS["minimap2"],
        TOOLS["samtools"],
        TOOLS["seqtk"],
        TOOLS["R"]
    params:
        ref=config[gencode_version]['ref'],
        rscript=join(SCRIPTSDIR,"pull_gff_classification.R"),
        gene_id=lambda wildcards: f"{config[wildcards.gene]['id']}", # Gene ID
    threads: 
        getthreads("pull_genes")
    shell:
        """
        # Get pigeon classifications and GFF for gene of interest / sample
        Rscript {params.rscript} {params.gene_id} {input.classification_filtered} {input.gff} {output.classification} {output.gff}

        while read line
            do 
                grep $line'|' {input.fasta} >> {output.tmp_list}
            done < <(cut -f1 {output.classification} | grep -v "isoform")
        sed -e 's/^>//' {output.tmp_list} > {output.id_list}

        # Filter fasta based on isoform ID
        seqtk subseq {input.fasta} {output.id_list} > {output.fasta}

        # Map fasta to genome using minimap2
        minimap2 -ax splice -t {threads} -uf --secondary=no -C5 {params.ref} {output.fasta} | samtools view -b -o {output.bam}

        # Sort and index bam
        samtools sort -o {output.sorted_bam} -@ {threads} {output.bam}
        samtools index {output.sorted_bam}
        """