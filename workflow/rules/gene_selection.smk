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
	    ref=config[{gencode_version}]['ref'],
		rscript=join(SCRIPTSDIR,"pull_gff_classification.R"),
		sample="{sample}",
		gene="{gene}"
    threads: 
        getthreads("pull_genes")
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
        minimap2 -ax splice -t {threads} -uf --secondary=no -C5 {params.ref} genes_of_interest/$SAMPLE\.$GENE\.fasta | samtools view -b -o {output.bam}

        # Sort and index bam
        samtools sort -o {output.sorted_bam} -@ {threads} {output.bam}
        samtools index {output.sorted_bam}
		"""

rule convert_majiq:
input:
output:
params:
shell:
# Convert GFF to GTF using gffread
for file in *.sorted.filtered_lite.gff; do /data/khanlab3/gb_omics/software/gffread/gffread -E $file -T -o $(basename $file .gff).gtf; done

convert_gtf_majiq.R