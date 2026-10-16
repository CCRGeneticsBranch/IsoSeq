# IsoSeq Pipeline

Snakemake pipeline for PacBio IsoSeq long-read transcriptome preprocessing and downstream analysis. Following the [PacBio Iso-Seq documentation](https://isoseq.how/), the pipeline removes primers and refines full-length reads, clusters and collapses transcripts against the genome, and classifies and filters isoforms with Pigeon. It also calls fusions with PBFusion and LongGF, maps EWSR1 fusion isoforms to custom fusion genomes, extracts isoforms for genes of interest, and converts the results to a MAJIQ-compatible format.

## Pipeline steps

| Step | Rule | Tool | Description |
|---|---|---|---|
| 1 | `link_ccs` | | Link CCS BAMs (and `.pbi` indexes) into `results/ccs_reads/` with sample names |
| 2 | `cluster` | lima, IsoSeq3 `refine`, `cluster` | Remove primers, keep full-length non-chimeric (FLNC) reads with poly(A) tails, and cluster into high-quality isoforms |
| 3 | `collapse` | pbmm2, IsoSeq3 `collapse` | Map clustered isoforms to the genome and collapse redundant isoforms by exon structure |
| 4 | `sort_annotation` | Pigeon `prepare` | Sort and index the GENCODE annotation |
| 5 | `pigeon` | Pigeon `prepare`, `classify`, `filter`, `report` | Classify isoforms against the annotation, filter artifacts, and report gene saturation |
| 6 | `pbfusion` | PBFusion | Detect fusion isoforms with optimized parameters |
| 7 | `longgf` | samtools, LongGF | Detect fusion isoforms |
| 8 | `pull_fusion`, `map_fusion` | seqtk, minimap2, samtools | Extract EWSR1 fusion isoforms from PBFusion calls and map them to each custom fusion genome |
| 9 | `pull_genes` | R, seqtk, minimap2, samtools | Extract classification, structure, sequence, and alignment for each gene of interest |
| 10 | `convert_majiq` | gffread, R | Convert filtered isoforms to GTF with MAJIQ-compatible attributes and per-isoform counts |

## Requirements

- Runs on NIH [Biowulf](https://hpc.nih.gov/) (SLURM and environment modules).
- Must have access to dependencies listed in [config/config.yaml](config/config.yaml).
- Requires the local `pacbio` conda environment, activated through `myconda` on pipeline initiation. The environment is not included in this repository; for setup, see the [Biowulf Conda instructions](https://hpc.nih.gov/apps/conda.html). It requires the Bioconda packages [isoseq3](https://bioconda.github.io/recipes/isoseq3/README.html), [lima](https://anaconda.org/bioconda/lima), [pbmm2](https://anaconda.org/bioconda/pbmm2), [pbpigeon](https://anaconda.org/bioconda/pbpigeon), [pbfusion](https://anaconda.org/bioconda/pbfusion), [LongGF](https://anaconda.org/bioconda/longgf), and [gffread](https://anaconda.org/bioconda/gffread).

## Installation

```bash
git clone https://github.com/CCRGeneticsBranch/IsoSeq.git
```

## Quick start

The pipeline is designed to be run in a separate directory (/path/to/my_run).

```bash
# 1. Create a work folder with template config files
bash IsoSeq/run_isoseq --runmode=init --workdir=/path/to/my_run [--samplesheet=/path/to/samples.tsv]

# 2. Edit /path/to/my_run/config/config.yaml and config/samples.tsv (see Inputs)

# 3. Dry run (writes /path/to/my_run/dryrun.log)
bash IsoSeq/run_isoseq --runmode=dryrun --workdir=/path/to/my_run

# 4. Submit to SLURM
bash IsoSeq/run_isoseq --runmode=run --workdir=/path/to/my_run

# Unlock the work folder after an interrupted run
bash IsoSeq/run_isoseq --runmode=unlock --workdir=/path/to/my_run
```

`-h` shows all options; `-v` prints the pipeline and launcher versions. Per-rule SLURM resources are in `WORKDIR/config/cluster.yaml` (copied from [resources/cluster.yaml](resources/cluster.yaml) at `init`).

## Inputs

### Sample sheet (`config/samples.tsv`)

| Column | Description | Example |
|---|---|---|
| `Sample` | Sample name | `sample1` |
| `CCS` | Absolute path to the PacBio CCS (HiFi) BAM; a `.pbi` index must exist alongside it | `/path/to/sample1.ccs.bam` |

### Configuration (`config/config.yaml`)

| Key | Default | Description (source) | Edit? |
|---|---|---|---|
| `gencode_version` | `v28` | GENCODE release used for the reference genome and annotation (`v28`, `v41`; [GENCODE human releases](https://www.gencodegenes.org/human/releases.html)) | As needed |
| `v28`, `v41`: `ref`, `gtf` | `/data/khanlab3/gb_omics/resources/` | GENCODE GRCh38 primary assembly FASTA and annotation GTF | No |
| `primers` | `PIPELINE_HOME/resources/primers.fasta` | [IsoSeq primer sequences](https://isoseq.how/clustering/cli-workflow.html) | No |
| `fusion_genome` | `[ewsr1_fli1, ewsr1_erg]` | Custom fusion genomes for mapping EWSR1 fusion isoforms | As needed |
| `<fusion_genome>: fusion_ref` | `PIPELINE_HOME/resources/fusion_genome/*.fa` | Custom fusion genome FASTA | No |
| `genes_of_interest` | `[EWSR1, FLI1]` | Gene names to extract; each must have an ID under the selected GENCODE release | As needed |
| `v28`, `v41`: `genes_of_interest.<gene>.id` | Versioned Ensembl gene IDs | Gene IDs associated with each GENCODE release; `gencode_version` selects the matching IDs | No |
| `workdir`, `manifest`, `tools`, `cluster` | `WORKDIR/...`, `PIPELINE_HOME/...` | Work folder and launcher/repository paths (filled by `init`) | No |

**Note:** The custom fusion genomes consist of the full-length fusion partner gene sequences concatenated together in the correct orientation (same coordinates in Gencode v28 and v41). 

## Outputs

All paths are relative to the work folder.

| Path | Description |
|---|---|
| `results/cluster/<sample>.flnc.bam`, `.clustered.hq.bam`, `.clustered.hq.fasta.gz` | FLNC reads and high-quality clustered isoforms |
| `results/collapse/<sample>.hq.mapped.bam`, `.hq.collapsed.gff`, `.fasta`, `.flnc_count.txt` | Mapped isoforms, collapsed isoform structures, sequences, and FLNC counts |
| `results/classify/<sample>_classification.txt`, `_classification.filtered_lite_classification.txt` | Pigeon isoform classifications (all and filtered) |
| `results/classify/<sample>.hq.collapsed.sorted.filtered_lite.gff`, `<sample>.saturation.txt` | Filtered isoforms and gene saturation report |
| `results/pbfusion/<sample>.breakpoints.groups.bed` | PBFusion fusion calls |
| `results/longgf/<sample>.log` | LongGF fusion calls |
| `results/fusion/<sample>_ewsr1.list`, `.fasta`, `<sample>.<fusion_genome>.sorted.bam` | EWSR1 fusion isoform IDs, sequences, and alignments to each fusion genome |
| `results/genes_of_interest/<sample>.<gene>.*` | Classification, GFF, FASTA, and sorted BAM for each gene of interest |
| `results/majiq/<sample>.sorted.filtered_lite.attributes.gtf`, `<sample>.counts.tsv` | MAJIQ-compatible GTF and per-isoform counts |
| `results/pipeline_version.txt` | Pipeline version used for the run |
| `runinfo.yaml`, `snakemake.log`, `logs/` | Run information, Snakemake log and per-rule SLURM logs |

## Tools and references

Modules are listed in [resources/tools.yaml](resources/tools.yaml); PacBio tools come from the `pacbio` conda environment.

| Tool | Module | Used in | Reference |
|---|---|---|---|
| Snakemake | `snakemake/7.19.1` | workflow | Mölder et al. 2021, *F1000Research* [doi:10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2) |
| lima | `pacbio` env | `cluster` | [PacificBiosciences/barcoding](https://github.com/PacificBiosciences/barcoding) |
| IsoSeq3 | `pacbio` env | `cluster`, `collapse` | [PacificBiosciences/IsoSeq](https://github.com/PacificBiosciences/IsoSeq) |
| pbmm2 | `pacbio` env | `collapse` | [PacificBiosciences/pbmm2](https://github.com/PacificBiosciences/pbmm2) |
| Pigeon | `pacbio` env | `sort_annotation`, `pigeon` | [PacificBiosciences/pigeon](https://github.com/PacificBiosciences/pigeon) |
| PBFusion | `pacbio` env | `pbfusion` | [PacificBiosciences/pbfusion](https://github.com/PacificBiosciences/pbfusion) |
| LongGF | `pacbio` env | `longgf` | Liu et al. 2020, *BMC Genomics* [doi:10.1186/s12864-020-07207-4](https://doi.org/10.1186/s12864-020-07207-4) |
| gffread | `pacbio` env | `convert_majiq` | Pertea & Pertea 2020, *F1000Research* [doi:10.12688/f1000research.23297.2](https://doi.org/10.12688/f1000research.23297.2) |
| minimap2 | `minimap2` | `map_fusion`, `pull_genes` | Li 2018, *Bioinformatics* [doi:10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191) |
| SAMtools | `samtools` | `longgf`, `map_fusion`, `pull_genes` | Danecek et al. 2021, *GigaScience* [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008) |
| seqtk | `seqtk` | `pull_fusion`, `pull_genes` | [lh3/seqtk](https://github.com/lh3/seqtk) |
| R | `R` | `pull_genes`, `convert_majiq` | [R Project](https://www.r-project.org/) |

## License

Original pipeline code and documentation are licensed under MIT; see [LICENSE](LICENSE). Third-party software, reference data and bundled sequences retain their own terms. The MIT license does not grant redistribution rights for external or bundled third-party resources.

## Versioning

The pipeline version is in [VERSION](VERSION) and changes are listed in [CHANGELOG.md](CHANGELOG.md). Each run records the pipeline version, git commit and launcher version in `runinfo.yaml` and `results/pipeline_version.txt`. The launch script is a copy of [biowulf-snakemake-launcher](https://github.com/CCRGeneticsBranch/biowulf-snakemake-launcher) (`bin/launcher.sh`).

## Contact

Erica Pehrsson, Genetics Branch, Center for Cancer Research, National Cancer Institute. Please report problems through GitHub issues.
