from contig_cds import ContigCDS
import itertools
from kmers import Kmers
import pandas as pd
import sys


class TrinucleotideFrequency(Kmers):
    def __post_init__(self):
        super().__post_init__()
        self.kmers = sorted(
            ["".join(p) for p in itertools.product(self.NUCLEOTIDES, repeat=self.k)]
        )
        self.frequency_table = pd.DataFrame(index=self.kmers)

    def count(self, sequence, counter):
        # Sliding window of 1
        for i in range(len(sequence) - self.k + 1):
            trimer = str(sequence[i : i + self.k])

            # Skip trimer with ambiguous bases
            if any(base not in "ATGC" for base in trimer):
                continue

            # Add count to seen trimer
            counter[trimer] += 1


if __name__ == "__main__":
    proteins_file = sys.argv[1]
    contigs_file = sys.argv[2]
    output_file = sys.argv[3]

    contig_cds = ContigCDS.parse_file(proteins_file)
    trinucleotides = TrinucleotideFrequency.set_k(3)
    trinucleotides.fill_frequency_table(contigs_file, contig_cds)
    trinucleotides.to_tsv(output_file)
