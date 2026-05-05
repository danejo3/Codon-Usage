from contig_cds import ContigCDS
import itertools
from kmers import Kmers
import pandas as pd
import sys


class CodonUsage(Kmers):
    def __post_init__(self):
        super().__post_init__()
        self.kmers = sorted(
            ["".join(p) for p in itertools.product(self.NUCLEOTIDES, repeat=self.k)]
        )
        self.frequency_table = pd.DataFrame(index=self.kmers)

    def count(self, sequence, counter):
        # Sliding window of 3
        for i in range(0, len(sequence) - self.k + 1, self.k):
            codon = sequence[i : i + self.k]

            # Skip codon with ambiguous bases
            if any(base not in "ATGC" for base in codon):
                continue

            # Add count to seen codon
            counter[codon] += 1


if __name__ == "__main__":
    proteins_file = sys.argv[1]
    contigs_file = sys.argv[2]
    output_file = sys.argv[3]

    contig_cds = ContigCDS.parse_file(proteins_file)
    codon_usage = CodonUsage.set_k(3)
    codon_usage.fill_frequency_table(contigs_file, contig_cds)
    codon_usage.to_tsv(output_file)
