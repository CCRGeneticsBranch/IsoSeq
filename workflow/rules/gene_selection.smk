rule pull_genes:
    """
    Pull genes of interest from IsoSeq results
    """
    input:
        fasta=join(RESULTSDIR,"collapse","{sample}.hq.collapsed.fasta") # Full fasta headers per isoform
    output:
        classification=join(RESULTSDIR,"genes_of_interest","{sample}.{gene}.classification.filtered_lite_classification.txt"),
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
        gene="{gene}", # Gene name
        gene_id=lambda wildcards: f"{config[wildcards.gene]['id']}", # Gene ID
        sample="{sample}"
    threads: 
        getthreads("pull_genes")
    shell:
        """
        # Get pigeon classifications and GFF for gene of interest / sample
        Rscript {params.rscript} {params.gene} {params.gene_id} {params.sample}

        while read line
            do 
                grep $line"|" {input.fasta} >> {output.tmp_list}
            done < <(cut -f1 {output.classification})
        sed -e 's/^>//' {output.tmp_list} > {output.id_list}

        # Filter fasta based on isoform ID
        seqtk subseq {input.fasta} {output.id_list} > {output.fasta}

        # Map fasta to genome using minimap2
        minimap2 -ax splice -t {threads} -uf --secondary=no -C5 {params.ref} {output.fasta} | samtools view -b -o {output.bam}

        # Sort and index bam
        samtools sort -o {output.sorted_bam} -@ {threads} {output.bam}
        samtools index {output.sorted_bam}
        """

rule convert_majiq:
    input:
        gff=join(RESULTSDIR,"classify","{sample}.sorted.filtered_lite.gff"),
        classification=join(RESULTSDIR,"classify","{sample}.classification.filtered_lite_classification.txt"),
    output:
        gtf=join(RESULTSDIR,"majiq","{sample}.sorted.filtered_lite.gtf"),
        gtf_attr=join(RESULTSDIR,"majiq","{sample}.sorted.filtered_lite.attributes.gtf"),
        counts=join(RESULTSDIR,"majiq","{sample}.counts.tsv")
    envmodules:
        TOOLS["R"]
    params:
        rscript=join(SCRIPTSDIR,"convert_gtf_majiq.R"),
        sample="{sample}"
    shell:
        """
        # Convert GFF to GTF using gffread
        /data/khanlab3/gb_omics/software/gffread/gffread -E {input.gff} -T -o {output.gtf}

        Rscript {params.rscript} {params.sample}
        """