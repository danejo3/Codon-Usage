from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
from dataclasses import dataclass, field
import numpy as np
import pandas as pd
from typing import ClassVar


@dataclass
class Kmers:
    k: int
    kmers: list[str] = field(init=False)
    frequency_table: np.ndarray = field(init=False)

    NUCLEOTIDES: ClassVar[list[str]] = ["A", "T", "G", "C"]

    def __post_init__(self):
        pass

    @classmethod
    def set_k(cls, k):
        return cls(k=k)

    def fill_frequency_table(self, contigs_file, contig_cds):
        all_counts = {}
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

            # Save all frequency counts for contig
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

    def get_cds(self, sequence, start, end, direction):
        if direction == 1:
            return sequence[start - 1 : end]
        else:
            return str(Seq(sequence[start - 1 : end]).reverse_complement())

    def to_tsv(self, output_file):
        self.frequency_table.to_csv(output_file, sep="\t")
