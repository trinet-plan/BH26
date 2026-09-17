"""Indexed FASTA access and VCF normalization without loading a genome into memory."""

from pathlib import Path

from acmg_pipeline.automated_core.models import Variant


class FastaReference:
    """Read an uncompressed FASTA with its existing samtools-compatible .fai index."""

    def __init__(self, path):
        self.path = Path(path)
        index_path = Path(str(self.path) + ".fai")
        self.index = {}
        for line in index_path.read_text(encoding="utf-8").splitlines():
            fields = line.split("\t")
            if len(fields) < 5:
                raise ValueError("Invalid FASTA index")
            name, length, offset, bases, width = fields[:5]
            values = tuple(map(int, (length, offset, bases, width)))
            if values[0] <= 0 or values[1] < 0 or values[2] <= 0 or values[3] < values[2]:
                raise ValueError("Invalid FASTA index dimensions")
            key = name.removeprefix("chr")
            if key in self.index:
                raise ValueError(f"Ambiguous reference sequence: {key}")
            self.index[key] = values

    def sequence(self, chrom, start, end):
        """Return 1-based inclusive reference interval, refusing out-of-bounds requests."""
        chrom = chrom.removeprefix("chr")
        if chrom not in self.index:
            raise ValueError(f"Chromosome not found in reference: {chrom}")
        length, offset, bases, width = self.index[chrom]
        if not 1 <= start <= end <= length:
            raise ValueError("Reference interval out of bounds")
        first = offset + ((start - 1) // bases) * width + (start - 1) % bases
        last = offset + ((end - 1) // bases) * width + (end - 1) % bases
        with self.path.open("rb") as stream:
            stream.seek(first)
            result = stream.read(last - first + 1).replace(b"\r", b"").replace(b"\n", b"")
        value = result.decode("ascii").upper()
        if len(value) != end - start + 1 or any(base not in "ACGT" for base in value):
            raise ValueError("Reference interval is unavailable or contains ambiguous bases")
        return value


def normalize(variant, reference):
    """Verify REF, left-align indels, then return minimal anchored VCF alleles."""
    pos, ref, alt = variant.pos, variant.ref, variant.alt
    observed = reference.sequence(variant.chrom, pos, pos + len(ref) - 1)
    if observed != ref:
        raise ValueError(f"REF_MISMATCH: expected {observed}, received {ref}")
    # Strip suffixes. If an allele empties, extend left so an indel can rotate through repeats.
    while ref[-1] == alt[-1]:
        if min(len(ref), len(alt)) == 1:
            if pos == 1:
                break
            base = reference.sequence(variant.chrom, pos - 1, pos - 1)
            ref, alt = base + ref, base + alt
            pos -= 1
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return Variant(variant.assembly, variant.chrom, pos, ref, alt)
