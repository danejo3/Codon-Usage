import gzip
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


SAMPLES = ["r1", "r2"][: len(config["config"]["reads"])]
SCRIPTS_DIR = config["config"]["scripts_dir"]
WEIGHTS = ["before_weights", "after_weights"]
FEATURES = [
    "tetranucleotide_frequency_on_contigs",
    "trinucleotide_frequency_on_cds",
    "codon_usage_on_cds",
    "rscu_on_cds",
]


rule all:
    input:
        "seq/map_reads_to_references/coverage_depth.tsv",
        "seq/map_reads_to_references/mapped.tsv",
        expand(
            "analysis/{weights}/{feature}/2d_pca.png",
            weights=WEIGHTS,
            feature=FEATURES,
        ),
        expand("analysis/{weights}/combined_pca_plots.png", weights=WEIGHTS),
        expand("analysis/{weights}/combined_silhouette_scores.txt", weights=WEIGHTS),


rule copy_input_reads:
    input:
        reads=config["config"]["reads"],
    output:
        reads=expand("seq/input/{sample}.fastq.gz", sample=SAMPLES),
    run:
        Path("seq").mkdir(parents=True, exist_ok=True)
        for i, f in enumerate(input.reads):
            out = output.reads[i]
            if f.endswith(".gz"):
                shell(f"ln -s {f} {out}")
            else:
                shell(f"gzip -c {f} > {out}")


rule fastp:
    input:
        reads=rules.copy_input_reads.output.reads,
    output:
        reads=expand("seq/fastp/{sample}.fastq.gz", sample=SAMPLES),
    threads: 8
    params:
        html_report="seq/fastp/fastp.html",
        json_report="seq/fastp/fastp.json",
    run:
        if len(input.reads) == 1:
            shell(
                "fastp -i {input.reads} -o {output.reads} -w {threads} --html {params.html_report} --json {params.json_report}"
            )
        else:
            shell(
                "fastp -i {input.reads[0]} -I {input.reads[1]} -o {output.reads[0]} -O {output.reads[1]} -w {threads} --html {params.html_report} --json {params.json_report}"
            )


rule combine_references:
    input:
        references=config["config"]["references"],
    output:
        references="seq/input/references.fasta",
    shell:
        """
        cat {input.references} > {output.references}
        """


rule map_reads_to_references:
    input:
        references=rules.combine_references.output.references,
        reads=rules.fastp.output.reads,
    output:
        aln="seq/map_reads_to_references/aln.bam",
        stats="seq/map_reads_to_references/stats.tsv",
        coverage_depth="seq/map_reads_to_references/coverage_depth.tsv",
        mapped="seq/map_reads_to_references/mapped.tsv",
    threads: 8
    shell:
        """
        minimap2 -ax sr {input.references} {input.reads} -t {threads} | samtools sort -@ {threads} -o {output.aln}
        samtools index {output.aln}
        samtools idxstats {output.aln} > {output.stats}
        samtools coverage {output.aln} > {output.coverage_depth}
        samtools view -F 256 {output.aln} | cut -f1,3 > {output.mapped}
        """


rule assembly:
    input:
        reads=rules.fastp.output.reads,
    output:
        contigs="analysis/assembly/contigs.fasta",
    threads: 8
    params:
        megahit_outdir="analysis/assembly/megahit",
        spades_outdir="analysis/assembly/spades",
    run:
        if len(input.reads) == 1:
            shell("megahit -r {input.reads} -o {params.megahit_outdir} -t 1")
            shell(
                "seqkit seq -m 2500 {params.megahit_outdir}/final.contigs.fa > {output.contigs}"
            )
        else:
            shell(
                "spades.py -1 {input.reads[0]} -2 {input.reads[1]} -o {params.spades_outdir} --meta -t {threads}"
            )
            shell(
                "seqkit seq -m 2500 {params.spades_outdir}/contigs.fasta > {output.contigs}"
            )


rule map_reads_to_contigs:
    input:
        contigs=rules.assembly.output.contigs,
        reads=rules.fastp.output.reads,
    output:
        aln="analysis/assembly/map_reads_to_contigs/aln.bam",
        stats="analysis/assembly/map_reads_to_contigs/stats.tsv",
        coverage_depth="analysis/assembly/map_reads_to_contigs/coverage_depth.tsv",
        mapped="analysis/assembly/map_reads_to_contigs/mapped.tsv",
    threads: 8
    shell:
        """
        minimap2 -ax sr {input.contigs} {input.reads} -t {threads} | samtools sort -@ {threads} -o {output.aln}
        samtools index {output.aln}
        samtools idxstats {output.aln} > {output.stats}
        samtools coverage {output.aln} > {output.coverage_depth}
        samtools view -F 256 {output.aln} | cut -f1,3 > {output.mapped}
        """


rule map_contigs_to_references:
    input:
        contigs=rules.assembly.output.contigs,
        references=rules.combine_references.output.references,
    output:
        aln="analysis/assembly/map_contigs_to_references/aln.bam",
        stats="analysis/assembly/map_contigs_to_references/stats.tsv",
        coverage_depth="analysis/assembly/map_contigs_to_references/coverage_depth.tsv",
        mapped="analysis/assembly/map_contigs_to_references/mapped.tsv",
    threads: 8
    shell:
        """
        minimap2 -ax asm5 {input.references} {input.contigs} -t {threads} | samtools sort -@ {threads} -o {output.aln}
        samtools index {output.aln}
        samtools idxstats {output.aln} > {output.stats}
        samtools coverage {output.aln} > {output.coverage_depth}
        samtools view -F 256 {output.aln} | cut -f1,3 > {output.mapped}
        """


rule contig_and_reference_colors:
    input:
        references=rules.combine_references.output.references,
        mapped=rules.map_contigs_to_references.output.mapped,
    output:
        contig_colors="analysis/assembly/map_contigs_to_references/contig_colors.tsv",
        reference_colors="analysis/assembly/map_contigs_to_references/reference_colors.tsv",
    run:
        species_map = {}
        species_counts = {}

        # Build species map
        with open(input.references) as f:
            for line in f:
                if line.startswith(">"):
                    parts = line[1:].strip().split()
                    accession = parts[0]
                    species = " ".join(parts[1:3])
                    species_map[accession] = species

                # Map accessions to species
        df = pd.read_csv(input.mapped, sep="\t", header=None)
        col = df.columns[1]
        df[col] = df[col].map(species_map).fillna("Unmapped")
        unique_species = sorted(df[col].unique())

        # Set color pallete
        color_palette = "tab10" if len(unique_species) <= 10 else "husl"
        palette = sns.color_palette(color_palette, len(unique_species))
        colors = {
            species: mcolors.to_hex(palette[i])
            for i, species in enumerate(unique_species)
        }
        # Force Unmapped values to black
        if "Unmapped" in unique_species:
            colors["Unmapped"] = "black"

            # Apply colors to dataframe
        ref_df = pd.DataFrame(list(colors.items()), columns=["species", "color"])
        ref_df.to_csv(output.reference_colors, sep="\t", index=False)
        df[col] = df[col].map(colors)
        df.to_csv(output.contig_colors, sep="\t", index=False)


rule gc_content_on_contigs:
    input:
        contigs=rules.assembly.output.contigs,
    output:
        gc_content="analysis/assembly/gc_content_on_contigs.tsv",
    shell:
        """
        seqkit fx2tab --name --only-id --gc {input.contigs} > {output.gc_content}
        """


rule tetranucleotides_on_contigs:
    input:
        contigs=rules.assembly.output.contigs,
    output:
        tetranucleotide_frequency="analysis/assembly/tetranucleotide_frequency_on_contigs.tsv",
    shell:
        """
        python {SCRIPTS_DIR}/tetranucleotide_frequency.py {input.contigs} {output.tetranucleotide_frequency}
        """


rule prodigal_on_contigs:
    input:
        contigs=rules.assembly.output.contigs,
    output:
        gbk="analysis/assembly/prodigal/contigs.gbk",
        proteins="analysis/assembly/prodigal/contigs.proteins.faa",
    shell:
        """
        prodigal -i {input.contigs} -o {output.gbk} -a {output.proteins}
        """


rule trinucleotide_on_cds:
    input:
        proteins=rules.prodigal_on_contigs.output.proteins,
        contigs=rules.assembly.output.contigs,
    output:
        trinucleotide_frequency="analysis/assembly/trinucleotide_frequency_on_cds.tsv",
    shell:
        """
        python {SCRIPTS_DIR}/trinucleotide_frequency.py {input.proteins} {input.contigs} {output.trinucleotide_frequency}
        """


rule codon_usage_on_cds:
    input:
        proteins=rules.prodigal_on_contigs.output.proteins,
        contigs=rules.assembly.output.contigs,
    output:
        codon_usage="analysis/assembly/codon_usage_on_cds.tsv",
    shell:
        """
        python {SCRIPTS_DIR}/codon_usage.py {input.proteins} {input.contigs} {output.codon_usage}
        """


rule rscu_on_cds:
    input:
        proteins=rules.prodigal_on_contigs.output.proteins,
        contigs=rules.assembly.output.contigs,
    output:
        rscu="analysis/assembly/rscu_on_cds.tsv",
    shell:
        """
        python {SCRIPTS_DIR}/rscu.py {input.proteins} {input.contigs} {output.rscu}
        """


def combine_tables(df, gc, cov, final):
    # Transpose
    kmer_df = df.T

    # Add GC content
    gc_data = gc.set_index(0)[1].to_dict()
    kmer_df["gc_content"] = kmer_df.index.map(gc_data)

    # Add coverage
    coverage_data = cov.set_index("#rname")["meandepth"].to_dict()
    kmer_df["coverage_depth"] = kmer_df.index.map(coverage_data)

    # Scale
    scaler_final = StandardScaler()
    scaled = scaler_final.fit_transform(kmer_df)
    scaled_df = pd.DataFrame(scaled, index=kmer_df.index, columns=kmer_df.columns)
    scaled_df.to_csv(final, sep="\t")


rule combine_feature_tables:
    input:
        feature=lambda wc: f"analysis/assembly/{wc.feature}.tsv",
        gc_content=rules.gc_content_on_contigs.output.gc_content,
        coverage_depth=rules.map_reads_to_contigs.output.coverage_depth,
    output:
        final="analysis/combined/{feature}.tsv",
    run:
        df = pd.read_csv(input.feature, sep="\t", index_col=0)
        gc = pd.read_csv(input.gc_content, sep="\t", header=None)
        cov = pd.read_csv(input.coverage_depth, sep="\t")
        combine_tables(df, gc, cov, output.final)


rule before_weights:
    input:
        combined=lambda wc: f"analysis/combined/{wc.feature}.tsv",
    output:
        final="analysis/before_weights/{feature}/final.tsv",
    run:
        df = pd.read_csv(input.combined, sep="\t", index_col=0)
        kmer_cols = df.columns[:-2]
        df[kmer_cols] = df[kmer_cols] / np.sqrt(len(kmer_cols))
        df.to_csv(output.final, sep="\t")


rule apply_weights:
    input:
        combined=lambda wc: f"analysis/combined/{wc.feature}.tsv",
    output:
        final="analysis/after_weights/{feature}/final.tsv",
    run:
        df = pd.read_csv(input.combined, sep="\t", index_col=0)
        kmer_cols = df.columns[:-2]
        df[kmer_cols] = df[kmer_cols] * 10 / np.sqrt(len(kmer_cols))
        df["gc_content"] *= 1
        df["coverage_depth"] *= 5
        df.to_csv(output.final, sep="\t")


rule plot_2d_pca:
    input:
        final="analysis/{weights}/{feature}/final.tsv",
        contig_colors=rules.contig_and_reference_colors.output.contig_colors,
        reference_colors=rules.contig_and_reference_colors.output.reference_colors,
    output:
        plot="analysis/{weights}/{feature}/2d_pca.png",
        score="analysis/{weights}/{feature}/score.txt",
    shell:
        """
        python {SCRIPTS_DIR}/2d_pca.py {input.final} {input.contig_colors} {input.reference_colors} {output.plot} "{wildcards.weights} - {wildcards.feature}" {output.score}
        """


rule combine_images:
    input:
        tetranucleotide="analysis/{weights}/tetranucleotide_frequency_on_contigs/2d_pca.png",
        trinucleotide="analysis/{weights}/trinucleotide_frequency_on_cds/2d_pca.png",
        codon="analysis/{weights}/codon_usage_on_cds/2d_pca.png",
        rscu="analysis/{weights}/rscu_on_cds/2d_pca.png",
    output:
        combined="analysis/{weights}/combined_pca_plots.png",
    shell:
        """
        magick {input} +append {output.combined}
        """


rule combine_silhouette_scores:
    input:
        tetranucleotide="analysis/{weights}/tetranucleotide_frequency_on_contigs/score.txt",
        trinucleotide="analysis/{weights}/trinucleotide_frequency_on_cds/score.txt",
        codon="analysis/{weights}/codon_usage_on_cds/score.txt",
        rscu="analysis/{weights}/rscu_on_cds/score.txt",
    output:
        combined="analysis/{weights}/combined_silhouette_scores.txt",
    shell:
        """
        echo -e "feature\tsilhouette_score" > {output.combined}
        echo -e "tetranucleotide\t$(cat {input.tetranucleotide})" >> {output.combined}
        echo -e "trinucleotide\t$(cat {input.trinucleotide})" >> {output.combined}
        echo -e "codon\t$(cat {input.codon})" >> {output.combined}
        echo -e "rscu\t$(cat {input.rscu})" >> {output.combined}
        """
