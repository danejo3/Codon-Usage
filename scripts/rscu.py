from Bio import SeqIO
from collections import Counter
from contig_cds import ContigCDS
import itertools
from kmers import Kmers
import pandas as pd
import sys


AA_TO_CODONS = {
    "F": ["TTT", "TTC"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "I": ["ATT", "ATC", "ATA"],
    "M": ["ATG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "Y": ["TAT", "TAC"],
    "H": ["CAT", "CAC"],
    "Q": ["CAA", "CAG"],
    "N": ["AAT", "AAC"],
    "K": ["AAA", "AAG"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "C": ["TGT", "TGC"],
    "W": ["TGG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
}


class RSCU(Kmers):
    def __post_init__(self):
        super().__post_init__()
        self.kmers = sorted(
            ["".join(p) for p in itertools.product(self.NUCLEOTIDES, repeat=self.k)]
        )
        self.frequency_table = pd.DataFrame(index=self.kmers)

    def fill_frequency_table(self, contigs_file, contig_cds):
        all_rscu = {}
        for record in SeqIO.parse(contigs_file, "fasta"):
            # Skip contig with no CDS
            if record.name not in contig_cds.contig_cds:
                continue

            # Loop through all CDS in contig and get frequency count
            counter = Counter()
            for cds in contig_cds.iter_cds(record.name):
                # Get CDS and count kmers
                sequence = self.get_cds(record.seq, cds.start, cds.end, cds.strand)
                self.count(sequence, counter)

            # Calculate RSCU for each codon
            rscu = self.compute_rscu(counter)

            # Save all RSCU for contig
            all_rscu[record.name] = rscu

        # Build DataFrame and normalize data to frequencies
        self.frequency_table = (
            pd.DataFrame.from_dict(all_rscu, orient="columns").fillna(0).astype("float64")
        )

    def count(self, sequence, counter):
        # Sliding window of 3
        for i in range(0, len(sequence) - self.k + 1, self.k):
            codon = sequence[i : i + self.k]

            # Skip codon with ambiguous bases
            if any(base not in "ATGC" for base in codon):
                continue

            # Add count to seen codon
            counter[codon] += 1

    def compute_rscu(self, counter):
        codon_counts = Counter({str(k): v for k, v in counter.items()})

        rscu = {}
        for aa, codons in AA_TO_CODONS.items():
            # Get counts for codons for this amino acid
            counts = [codon_counts.get(codon, 0) for codon in codons]
            total = sum(counts)
            n = len(codons)

            # Avoid division by zero
            if total == 0:
                for codon in codons:
                    rscu[codon] = 0.0
                continue

            # Apply RSCU formula
            expected = total / n
            for codon, count in zip(codons, counts):
                rscu[codon] = count / expected
        return rscu


if __name__ == "__main__":
    proteins_file = sys.argv[1]
    contigs_file = sys.argv[2]
    output_file = sys.argv[3]

    contig_cds = ContigCDS.parse_file(proteins_file)
    RSCU = RSCU.set_k(3)
    RSCU.fill_frequency_table(contigs_file, contig_cds)
    RSCU.to_tsv(output_file)
