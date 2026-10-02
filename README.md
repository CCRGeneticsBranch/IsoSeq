Runs the PacBio IsoSeq pipeline to:
1. Cluster, collapse, and classify transcripts
2. Call fusions using PBFusion and LongGF and map to a custom fusion genome
3. Extract annotation and sequence for genes of interest
4. Convert PacBio output to MAJIQ-compatible format

To run the pipeline, clone the repo from github. In an interactive session, initialize the analysis directory with:
```
  bash [path to repo]/IsoSeq/run_isoseq --runmode=init --workdir=[directory]
```
Edit the config/samples.tsv and config/config.yaml files, including GENCODE version, resource locations, and gene ENSEMBL IDs.
Perform a dryrun with:
```
  bash [path to repo]/IsoSeq/run_isoseq --runmode=dryrun --workdir=[directory]
```
If successful, launch the pipeline with:
```
  bash [path to repo]/IsoSeq/run_isoseq --runmode=run --workdir=[directory]
```
Note: must have environment "pacbio" installed per Biowulf instructions, including: isoseq3 lima pbmm2 pbpigeon pbfusion longgf gffread
