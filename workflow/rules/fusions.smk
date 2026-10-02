localrules: pbfusion, longgf, pull_fusion, map_fusion

rule pbfusion:
    """
    PacBio fusion detection
    """
    input:
        collapse=join(RESULTSDIR,"collapse","{sample}.hq.mapped.bam")
    output:
        breakpoints=join(RESULTSDIR,"pbfusion","{sample}.breakpoints.groups.bed")
    params:
        gtf=join(RESULTSDIR,f"gencode.{gencode_version}.annotation.sorted.gtf"),
        outdir=join(RESULTSDIR,"pbfusion","{sample}")
    threads: 
        getthreads("pbfusion")
    shell:
        """
        # Detect fusions
        pbfusion discover -c 5 -s 10000 -M 0 -t {threads} --gtf {params.gtf} --output-prefix {params.outdir} {input.collapse}
        """

rule longgf:
    """ 
    LongGF fusion detection
    """
    input:
        collapse=join(RESULTSDIR,"collapse","{sample}.hq.mapped.bam")
    output:
        sorted_bam=join(RESULTSDIR,"longgf","{sample}.hq.sorted.bam"),
        log=join(RESULTSDIR,"longgf","{sample}.log")
    envmodules:
        TOOLS["samtools"]
    params:
        gtf=config[gencode_version]['gtf']
    shell:
        """
        # Sort bam files by position
        samtools sort -n -o {output.sorted_bam} {input.collapse}

        # Detect fusions
        LongGF {output.sorted_bam} {params.gtf} 100 50 100 > {output.log}
        """

rule pull_fusion:
    """
	Get EWSR1 fusion transcript IDs and sequences
	"""
    input:
        cluster=join(RESULTSDIR,"cluster","{sample}.clustered.hq.fasta.gz"), # IsoSeq clustered HQ fasta
        breakpoints=join(RESULTSDIR,"pbfusion","{sample}.breakpoints.groups.bed")
    output:
        ids=join(RESULTSDIR,"fusion","{sample}_ewsr1.list"), # Fusion isoform IDs
        fasta=join(RESULTSDIR,"fusion","{sample}_ewsr1.fasta")
    envmodules:
        TOOLS["seqtk"]
    shell:
        """
        # Get list of EWSR1 fusion transcript IDs
        grep 'EWSR1' {input.breakpoints} | cut -f12 | cut -f1 -d';' | sed 's/\,/\\n/g' | sed 's/RN\=//g' > {output.ids}
        
        # Filter fasta file based on isoform ID
        seqtk subseq {input.cluster} {output.ids} > {output.fasta}
        """

rule map_fusion:
    """
    Map PBFusion fusion isoforms to custom genome
    """
    input:
        fasta=join(RESULTSDIR,"fusion","{sample}_ewsr1.fasta")
    output:
        bam=temp(join(RESULTSDIR,"fusion","{sample}.{fusion_genome}.bam")),
        sorted_bam=join(RESULTSDIR,"fusion","{sample}.{fusion_genome}.sorted.bam") # Fusion isoforms mapped to fusion genome
    envmodules:
        TOOLS["minimap2"],
        TOOLS["samtools"]
    params:
        ref=lambda wildcards: f"{config[wildcards.fusion_genome]['fusion_ref']}",
    shell:
        """
        # Map fasta to genome using minimap2
        minimap2 -ax splice -t {threads} -uf --secondary=no -C5 {params.ref} {input.fasta} | samtools view -b -o {output.bam}

        # Sort and index
        samtools sort -o {output.sorted_bam} {output.bam}
        samtools index {output.sorted_bam}
        """