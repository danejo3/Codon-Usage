# Codon-Usage

To separate metagenomic DNA sequences, many binning algorithms rely on genomic signals to cluster similar sequences that are derived from the same organism. Tetranucleotide frequency (TNF) has been a pivotal feature in metagenomic binning, as it captures genome-wide compositional patterns. In addition to TNF, codon usage bias (or optimization) represents another informative signal derived from protein-coding regions. This study presents a tool that enables users to explore how codon usage may provide additional discriminatory signal for sequence binning. The results demonstrate that codon usage frequency and relative synonymous codon usage, alongside GC content and coverage depth, both before and after feature weighting, produce consistent taxonomic clustering patterns using the HMP staggered dataset. While these findings have potential, several limitations remain, including the need for validation on a real metagenomic dataset across a diverse microbial community and an assessment of different clustering algorithms in non-linear feature spaces.

## Installation

Users can create a Conda environment containing all required dependencies to run the workflow. Alternatively, packages can be installed manually using the list provided in the `environment.yml` file.

```bash
# Download conda packages
conda create -n codon-usage -f environment.yml

# Download the Codon-Usage Github repo
git clone https://github.com/danejo3/Codon-Usage.git
cd Codon-Usage
```

## Usage

```
usage: cli_parser.py [-h] -r REFERENCES [REFERENCES ...] [-w WORKDIR] [-t THREADS] reads [reads ...]

A command-line interface for running Snakemake workflows that generate tetranucleotide and trinucleotide frequency profiles, codon usage statistics, and RSCU
plots.

positional arguments:
  reads                 input FASTQ files (max 2)

options:
  -h, --help            show this help message and exit
  -r REFERENCES [REFERENCES ...], --references REFERENCES [REFERENCES ...]
                        reference FASTA files
  -w WORKDIR, --workdir WORKDIR
                        output directory for workflow files (default: current directory)
  -t THREADS, --threads THREADS
                        number of threads to use (default: 1)
```

## Example of Running the Workflow

To run the workflow, provide the input FASTQ reads and reference FASTA files using the following example command:

```bash
python cli_parser.py sample/SRR172903.fastq.gz -w sandbox -t 8 -r references/Cereibacter_sphaeroides.fasta references/Escherichia_coli.fasta references/Methanobrevibacter_smithii.fasta references/Pseudomonas_aeruginosa.fasta references/Staphylococcus_aureus.fasta references/Staphylococcus_epidermidis.fasta references/Streptococcus_mutans.fasta
```

## Workflow Overview

The workflow automates key stages of the pipeline, including read processing, assembly, feature extraction, and downstream analysis. The Snakemake directed acyclic graph (DAG) below illustrates all rules included in the workflow.

<img width="620" height="461" alt="Snakemake workflow DAG" src="https://github.com/user-attachments/assets/7a9e062f-b30a-4e6f-acba-e89d7bdc7b12" />

## Example PCA Output Plots Before And After Feature Weighting

Principal component analysis (PCA) demonstrated that contig sequence composition was sufficient to separate contigs into genome representative clusters prior to feature weighting (Fig. 1A–D). Across tetranucleotide frequency (Fig. 1A), trinucleotide frequency (Fig. 1B), codon usage (Fig. 1C), and relative synonymous codon usage (RSCU; Fig. 1D), contigs formed either distinct linear or dense drop-like clustering patterns corresponding to their respective genomes. Notably, E. coli exhibited a long and strong linear distribution along PC1. In contrast, the remaining genomes produced more compact clusters with reduced dimensional dispersion; S. aureus, however, showed moderate dispersion along PC2.

<img width="6550" height="1500" alt="image" src="https://github.com/user-attachments/assets/5899aa33-f1ce-49d3-88d7-5b62d9432c62" />

Quantitative cluster evaluation further demonstrated that all feature sets had moderate discriminatory power prior to weighting. Unweighted tetranucleotide frequencies produced the highest silhouette score (0.441), followed closely by trinucleotide frequencies (0.426), RSCU (0.419), and codon usage (0.413). The relatively small differences among these values indicate that no single compositional representation outperformed the others. Instead, the results suggest that a broad composition in the genomes contains most of the taxonomic information needed.

<img width="6550" height="1500" alt="image" src="https://github.com/user-attachments/assets/d00e533d-733b-491d-a444-297f0665ed8b" />

Applying feature weights substantially improved cluster densities. PCA revealed strong taxonomic separation across all weighted feature sets (Fig. 2A–D). For tetranucleotide frequency (Fig. 2A), trinucleotide frequency (Fig. 2B), codon usage (Fig. 2C), and relative synonymous codon usage (RSCU; Fig. 2D), adding weights tightened the spatial distribution of contigs and shifted previously linear patterns into compact drop-like clusters. In particular, the E. coli cluster exhibited less spread along PC1 and adopted a more compact distribution. Comparable reductions in linear trajectories were observed for the other species as well.

Although weighted ordinations appeared more compact, quantitative clustering performance declined across all feature sets following weighting. Silhouette scores decreased for tetranucleotide frequencies (0.390), trinucleotide frequencies (0.373), codon usage (0.374), and RSCU (0.400), indicating that feature weighting reduced overall cluster separability despite improving apparent visual compactness. This discrepancy highlights an important distinction between visual cluster cohesion and formal clustering quality. The tighter distributions observed after weighting likely reflects the compression of dominant variance axes rather than increased separation among taxa. Subsequently, weighting may have reduced biologically informative signals while simultaneously emphasizing lower-variance features that contributed less effectively to taxonomic separation.

## Test Dataset and References

Users can recreate the plots shown above by downloading and using the HMP mock staggered dataset and their references here. https://osf.io/7sbtx/overview.
