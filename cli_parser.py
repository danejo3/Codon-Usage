from argparse import ArgumentParser
from pathlib import Path
import subprocess


DESCRIPTION = (
    "A command-line interface for running Snakemake workflows that generate "
    "tetranucleotide and trinucleotide frequency profiles, codon usage statistics, "
    "and RSCU plots."
)


parser = ArgumentParser(description=DESCRIPTION)
parser.add_argument(
    "reads",
    type=Path,
    nargs="+",
    help="input FASTQ files (max 2)",
)
parser.add_argument(
    "-r",
    "--references",
    type=Path,
    nargs="+",
    required=True,
    help="reference FASTA files",
)
parser.add_argument(
    "-w",
    "--workdir",
    type=Path,
    default=Path("."),
    help="output directory for workflow files (default: current directory)",
)
parser.add_argument(
    "-t",
    "--threads",
    type=int,
    default=1,
    help="number of threads to use (default: 1)",
)

args = parser.parse_args()
if len(args.reads) > 2:
    parser.error("You can provide at most 2 FASTQ files.")
config = {
    "reads": [str(f.resolve()) for f in args.reads],
    "references": [str(f.resolve()) for f in args.references],
    "scripts_dir": str(Path("scripts").resolve()),
}
command = [
    "snakemake",
    "-s",
    "Workflow.smk",
    "--directory",
    str(args.workdir.resolve()),
    "-c",
    str(args.threads),
    "--config",
    f"config={config}",
]

try:
    subprocess.run(command)
except subprocess.CalledProcessError as e:
    print("Snakemake failed with error:", e.stderr)
