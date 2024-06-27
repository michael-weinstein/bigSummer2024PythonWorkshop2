import dataclasses
import typing
import functools
import collections
def returnZero() -> int:
    return 0

@dataclasses.dataclass
class DNA:
    sequence: str

    def __post_init__(self):
        self.sequence = self.sequence.upper()

    @property
    def length(self) -> int:
        return len(self.sequence)

    def __str__(self) -> str:
        return self.sequence

    def __len__(self) -> int:
        return self.length

    @functools.cached_property
    def baseCounts(self) -> typing.Dict[str, int]:
        counter = collections.defaultdict(returnZero)
        for base in self.sequence:
            counter[base] += 1
        return counter

    @functools.cached_property
    def gcContent(self) -> float:
        gc = self.baseCounts["G"] + self.baseCounts["C"]
        return gc / self.length

@dataclasses.dataclass
class QualityString:
    qualityString: str

    def __post_init__(self):
        self.qualityString = self.qualityString.upper()

    @property
    def length(self) -> int:
        return len(self.qualityString)

    def __str__(self) -> str:
        return self.qualityString

    def __len__(self) -> int:
        return self.length

    @functools.cached_property
    def phredValues(self) -> typing.List[int]:
        return [ord(character) - 33 for character in self.qualityString]


@dataclasses.dataclass
class SamLine:
    qname: str
    flag: int
    rname: str
    pos: int
    mapq: int
    cigar: str
    rnext: str
    pnext: int
    tlen: int
    seq: DNA
    qual: QualityString


def readSamFile(filePath :str, dataLimit :int =0) -> typing.List[SamLine]:
    samFile = open(filePath, "r")
    samLines = []
    for line in samFile:
        if dataLimit and len(samLines) >= dataLimit:
            break
        if line.startswith("@"):
            continue
        line = line.strip()
        fields = line.split()
        fixedValues = fields[:11]
        tags = fields[11:]
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = fixedValues
        pos = int(pos)
        mapq = int(mapq)
        pnext = int(pnext)
        tlen = int(tlen)
        seq = DNA(seq)
        qual = QualityString(qual)
        samLines.append(SamLine(qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual))
    samFile.close()
    return samLines