from Bio import SeqIO
from collections import Counter
from dataclasses import dataclass, field
import itertools
from kmers import Kmers
import pandas as pd
import sys
from typing import Dict


@dataclass
class TetranucleotideFrequency(Kmers):
    canonical_map: Dict[str, str] = field(init=False)

    def __post_init__(self):
        super().__post_init__()
        tetramers = [
            "".join(p) for p in itertools.product(self.NUCLEOTIDES, repeat=self.k)
        ]

        # Generate tetramer canonical mapping
        canonical_map = {}
        for tetramer in tetramers:
            rc = self.reverse_complement(tetramer)
            canonical = min(tetramer, rc)
            canonical_map[tetramer] = canonical
        self.canonical_map = canonical_map

        # 136 canonical tetramers
        self.kmers = sorted(set(canonical_map.values()))
        self.frequency_table = pd.DataFrame(index=self.kmers)

    def reverse_complement(self, sequence):
        complement = str.maketrans("ATGC", "TACG")
        return sequence.translate(complement)[::-1]

    def fill_frequency_table(self, contigs_file):
        all_counts = {}
        for record in SeqIO.parse(contigs_file, "fasta"):
            # Get and save all frequency counts from contig
            counter = Counter()
            self.count(record.seq, counter)
            all_counts[record.name] = counter

        # Build DataFrame and normalize data to frequencies
        self.frequency_table = (
            pd.DataFrame.from_dict(all_counts, orient="columns")
            .fillna(0)
            .astype("int64")
        )
        self.frequency_table = self.frequency_table.div(
            self.frequency_table.sum(axis=0), axis=1
        )

    def count(self, sequence, counter):
        # Sliding window of 1
        for i in range(len(sequence) - self.k + 1):
            tetramer = str(sequence[i : i + self.k])

            # Skip tetramer with ambiguous bases
            if any(base not in "ATGC" for base in tetramer):
                continue

            # Add count to seen tetramer
            canonical = self.canonical_map[tetramer]
            counter[canonical] += 1


if __name__ == "__main__":
    contigs_file = sys.argv[1]
    output_file = sys.argv[2]

    tetranucleotides = TetranucleotideFrequency.set_k(4)
    tetranucleotides.fill_frequency_table(contigs_file)
    tetranucleotides.to_tsv(output_file)
