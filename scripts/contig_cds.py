from Bio import SeqIO
from dataclasses import dataclass
from collections import namedtuple, defaultdict
from typing import Dict


CDS = namedtuple("CDS", ["start", "end", "strand"])


@dataclass
class ContigCDS:
    contig_cds: Dict[str, list[CDS]]

    @classmethod
    def parse_file(cls, proteins_file):
        contig_cds = defaultdict(list)
        for record in SeqIO.parse(proteins_file, "fasta"):
            contig_id = record.id.rsplit("_", 1)[0]
            parts = list(map(str.strip, record.description.split("#")))
            _, start, end, strand, *_ = parts
            contig_cds[contig_id].append(CDS(int(start), int(end), int(strand)))
        return cls(contig_cds=contig_cds)

    def iter_cds(self, key):
        return iter(self.contig_cds.get(key, []))
