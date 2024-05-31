# Get list of samples
SAMPLE=$1

# Other lists
{fusion_genome}: ewsr1_fli1, ewsr1_erg

# Set directories
DATADIR=
RESULTSDIR=isoseq_pipeline
SCRIPTSDIR=/data/khanlab3/gb_omics/projects/caplen/ccrgb_13/scripts
version="v28"

# Make output directories
ccs_reads
if [ ! -d cluster ];then mkdir cluster ;fi
if [ ! -d collapse ];then mkdir collapse ;fi
if [ ! -d classify ];then mkdir classify ;fi
if [ ! -d pbfusion ];then mkdir pbfusion ;fi
if [ ! -d longgf ];then mkdir longgf ;fi
genes_of_interest
fusion/bam

# Load environment
source myconda
mamba activate isoseq

# Load references (config)
primers: primers.fasta
ref: /data/khanlab3/gb_omics/Erica/software/genomes/gencode_human_28/GRCh38.primary_assembly.genome.fa
gtf: /data/khanlab3/gb_omics/Erica/software/gencode/gencode.v28.annotation.gtf # Use {version}
fusion_ref: ../fusion_genome/ewsr1_fli1.fa OR ewsr1_erg

# Link samples from raw data directory with updated names

# Rules

rule sort_annotation:
    """ 
    Sort and index genome annotation file
    """
    output:
	    sorted_gtf=join(RESULTSDIR,"gencode.{version}.annotation.sorted.gtf")
    params:
	    gtf=config['gtf']
    shell:
        """
		pigeon sort {params.gtf} -o {output.sorted_gtf}
		pigeon index {output.sorted_gtf}
	    """

rule cluster:
    """
    IsoSeq3 clustering
    """
    input:
        ccs=join(RESULTSDIR,"ccs_reads","{sample}.ccs.bam")
    output:
        fl_primer=temp(join(RESULTSDIR,"cluster","{sample}.fl.NEB_5p--NEB_Clontech_3p.bam")),
		fl_primer_index=temp(join(RESULTSDIR,"cluster","{sample}.fl.NEB_5p--NEB_Clontech_3p.bam.pbi")),
		fl=join(RESULTSDIR,"cluster","{sample}.fl.bam"),
	    flnc=join(RESULTSDIR,"cluster","{sample}.flnc.bam"),
	    clustered=join(RESULTSDIR,"cluster","{sample}.clustered.bam")
    params:
	    threads=$((SLURM_CPUS_PER_TASK-4)),
		primers=config['primers']
    shell:
        """
	    # Lima for primer removal
	    lima {input.ccs} {params.primers} {output.fl} --isoseq --peek-guess -j {params.threads}
	    
	    # IsoSeq3 refine reads
	    isoseq3 refine {output.fl_primer} {params.primers} {output.flnc} --require-polya -j {params.threads}
	    	    
	    # IsoSeq3 cluster
	    isoseq3 cluster {output.flnc} {output.clustered} --verbose --use-qvs -j {params.threads}
		"""

rule collapse:
    """ 
    IsoSeq3 mapping and collapse
    """
    input:
	    clustered_hq=join(RESULTSDIR,"cluster","{sample}.clustered.hq.bam")
    output:
		collapse=join(RESULTSDIR,"collapse","{sample}.hq.mapped.bam"),
        gff=join(RESULTSDIR,"collapse","{sample}.hq.collapsed.gff")
    params:
	    ref=config['ref'],
		threads=$SLURM_CPUS_PER_TASK
    shell:
        """
        # Map transcripts to reference genome
        pbmm2 align --preset ISOSEQ --sort {input.clustered_hq} {params.ref} {output.collapse} -j {params.threads}
    
        # Collapse redundant transcripts into unique isoforms via exon structure
        isoseq3 collapse --do-not-collapse-extra-5exons {output.collapse} {output.gff} -j {params.threads}
	    """

rule pigeon:
    """
    Pigeon classification and filtering
    """
    input:
	    gff=join(RESULTSDIR,"collapse","{sample}.hq.collapsed.gff"),
		abundance=join(RESULTSDIR,"collapse","{sample}.hq.collapsed.abundance.txt")
    output:
	    gff_sorted=join(RESULTSDIR,"classify","{sample}.sorted.gff"),
		classification=join(RESULTSDIR,"classify","{sample}_classification.txt"),
		classification_filtered=join(RESULTSDIR,"classify","{sample}_classification.filtered_lite_classification.txt"),
		saturation=join(RESULTSDIR,"classify","{sample}.saturation.txt")
    params:
	    ref=config['ref'],
	    gtf=join(RESULTSDIR,"gencode.{version}.annotation.sorted.gtf"),
		threads=$SLURM_CPUS_PER_TASK
		sample="{sample}"
    shell:
        """
        # Sort transcript GFF
        pigeon sort {input.gff} -o {output.gff_sorted}
        
        # Classify isoforms
        pigeon classify {output.gff_sorted} {params.gtf} {params.ref} -j {params.threads} -d "classify" -o {params.sample} \
            --gene-id --fl {input.abundance}

        # Filter isoforms
        pigeon filter {output.classification} --isoforms {output.gff_sorted} -j {params.threads}

        # Report gene saturation
        pigeon report {output.classification_filtered} {output.saturation} -j {params.threads}
        """

rule pbfusion:
    """
    PacBio fusion detection
    """
    input:
	    collapse=join(RESULTSDIR,"collapse","{sample}.hq.mapped.bam")
    output:
	    breakpoints=join(RESULTSDIR,"fusion","{sample}.breakpoints.groups.bed")
    params:
	    ref=config['ref'],
	    gtf=join(RESULTSDIR,"gencode.{version}.annotation.sorted.gtf"),
		threads=$SLURM_CPUS_PER_TASK,
		sample=join(RESULTSDIR,"fusion","{sample}")
    shell:
	    """
        # Detect fusions
        pbfusion discover -c 5 -s 10000 -M 0 -t {params.threads} --gtf {params.gtf} --output-prefix {params.sample} {input.collapse}
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
    params:
	    gtf=config['gtf']
    shell:
        """
        # Sort bam files by position
        samtools sort -n -o {output.sorted_bam} {input.collapse}

        # Detect fusions
        LongGF {output.sorted_bam} {params.gtf} 100 50 100 > {output.log}
        """

rule map_fusion:
    """
    Map fusion isoforms to custom genome
	"""
    input:
	    cluster=join(RESULTSDIR,"cluster","{sample}.clustered.hq.fasta.gz"), # IsoSeq clustered HQ fasta
		breakpoints=join(RESULTSDIR,"fusion","{sample}.breakpoints.groups.bed") # PBFusion output
    output:
	    ids=join(RESULTSDIR,"fusion","bam","{sample}_ewsr1.list") # Fusion isoform IDs
		fasta=join(RESULTSDIR,"fusion","bam","{sample}.fasta"),
		bam=temp(join(RESULTSDIR,"fusion","bam","{sample}_{fusion_genome}.bam")),
		sorted_bam=join(RESULTSDIR,"fusion","bam","{sample}_{fusion_genome}.bam") # Fusion isoforms mapped to fusion genome
    envmodules:
	    TOOLS["minimap2"],
		TOOLS["samtools"],
	    TOOLS["seqtk"]
    params:
	    ref=config['fusion_ref'], # Custom fusion reference
		threads=$SLURM_CPUS_PER_TASK
    shell:
	    """
		# Get list of EWSR1 fusion transcript IDs
		grep 'EWSR1' {input.breakpoints} | cut -f12 | cut -f1 -d";" | sed 's/\,/\n/g' | sed 's/RN\=//g' > {output.ids}
		
		# Filter fasta file based on isoform ID
		seqtk subseq {input.cluster} {output.ids} > {output.fasta}

		# Map fasta to genome using minimap2
		minimap2 -ax splice -t {params.threads} -uf --secondary=no -C5 {params.ref} {output.fasta} | samtools view -b -o {output.bam}

		# Sort and index
		samtools sort -o {output.sorted_bam} {output.bam}
		samtools index {output.sorted_bam}
		"""

rule pull_genes:
    """
    Usage: bash get_gene_isoform.sh [Gene name of interest] [Gene ID of interest] [Sample]
    
    EWSR1 - "ENSG00000182944.17"
    FLI1 - "ENSG00000151702.16"
    CD99 - "ENSG00000002586.19"
    PURPL - "ENSG00000250337.5"
	"""
    input:
	    fasta=join(RESULTSDIR,"collapse","{sample}.hq.collapsed.fasta") # Full fasta headers per isoform
	    GENE=$1 # Gene of interest
	    ID=$2
    output:
	    bam=temp(join(RESULTSDIR,"genes_of_interest","{sample}.{gene}.bam"))
		sorted_bam=join(RESULTSDIR,"genes_of_interest","{sample}.{gene}.sorted.bam")
    envmodules:
	    TOOLS["minimap2"],
		TOOLS["samtools"],
		TOOLS["seqtk"],
		TOOLS["R"]
    params:
	    ref=config['ref'],
		rscript=join(SCRIPTSDIR,"pull_gff_classification.R"),
		threads=$SLURM_CPUS_PER_TASK,
		sample="{sample}",
		gene="{gene}"
    shell:
        """
		# Get pigeon classifications and GFF for gene of interest / sample
        Rscript {params.rscript} {params.gene} $ID {params.sample}

        while read line
        	do 
		        grep $line"|" {input.fasta} >> genes_of_interest/$SAMPLE\.$GENE\.tmp
	        done < <(cut -f1 genes_of_interest/$SAMPLE\.$GENE\_classification.filtered_lite_classification.txt)
        sed -e 's/^>//' genes_of_interest/$SAMPLE\.$GENE\.tmp > genes_of_interest/$SAMPLE\.$GENE\.list
        rm genes_of_interest/$SAMPLE\.$GENE\.tmp

        # Filter fasta based on isoform ID
        seqtk subseq {input.fasta} genes_of_interest/$SAMPLE\.$GENE\.list > genes_of_interest/$SAMPLE\.$GENE\.fasta
        rm genes_of_interest/$SAMPLE\.$GENE\.list

            # Map fasta to genome using minimap2
        minimap2 -ax splice -t {params.threads} -uf --secondary=no -C5 {params.ref} genes_of_interest/$SAMPLE\.$GENE\.fasta | samtools view -b -o {output.bam}

        # Sort and index bam
        samtools sort -o {output.sorted_bam} -@ {params.threads} {output.bam}
        samtools index {output.sorted_bam}
		"""

rule convert_majiq:
input:
output:
params:
shell:
# Convert GFF to GTF using gffread (see OneNote)

convert_gtf_majiq.R