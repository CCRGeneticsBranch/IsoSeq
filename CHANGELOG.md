# Changelog

## 1.0.1
- Select gene-of-interest Ensembl IDs from the configured GENCODE release so gene extraction uses IDs matching the selected reference.

## 1.0.0
- Launch script replaced by the shared biowulf-snakemake-launcher (`bin/launcher.sh`); run modes: init, dryrun, unlock, run.
- Standardized repository layout (VERSION, CHANGELOG.md, .gitignore, resources/cluster.yaml, resources/tools.yaml).
- Pipeline version printed at start and written to results/pipeline_version.txt.
