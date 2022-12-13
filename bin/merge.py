import argparse
import logging
import sys
import tempfile
from pathlib import Path

from pybedtools import BedTool

logger = logging.getLogger()


class Coordinate:
    def __init__(self, chrom, pos, end):
        self._chrom = chrom
        self._pos = pos
        self._end = end

    @property
    def chrom(self):
        return self._chrom

    @property
    def pos(self):
        return self._pos

    @property
    def end(self):
        return self._end


class Variant:
    def __init__(self, chrom, pos, identifier, ref, alt, quality, filter, info, format_value, sample,
                 source=''):
        self._chrom = chrom
        self._pos = pos
        self._id = identifier
        self._ref = ref
        self._alt = alt
        self._quality = quality
        self._filter = filter
        self._info = info
        self._format = format_value
        self._sample = sample
        self._source = source

        self._set_type_len_end()
        self._set_confidence_and_sequence()

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, Variant):
            return self._chrom == other._chrom and self._pos == other._pos
        return False

    def _set_type_len_end(self):

        self._type = None
        self._len = 0
        self._end = None

        sv_len = len(self._alt) - len(self._ref)

        if self._info == 'INDEL':
            self._type = 'INDEL'

        if len(self._alt) == 1 and len(self._ref) == 1:
            self._type = 'SNP'
            self._calculated_type = 'SNP'
            self._end = self._pos

            self._info += f';CTYPE=SNP'

        info_list = self._info.split(';')

        for info in info_list:
            if 'SVTYPE=' in info:
                self._type = info.replace('SVTYPE=', '')
                continue

            if 'SVLEN=' in info:
                self._len = int(info.replace('SVLEN=', ''))
                continue

            if 'END=' in info:
                self._end = int(info.replace('END=', ''))
                continue

        if self.indel:
            self._calculated_type = 'INS' if sv_len > 0 else ('DEL' if sv_len < 0 else '')
            self._calculated_len = len(self._alt) - len(self._ref)
            self._calculated_end = self._pos + abs(sv_len)

            self._info += f';CTYPE={self._calculated_type};CLEN={self._calculated_len};CEND={self._calculated_end}'

    def _set_confidence_and_sequence(self):
        if self.indel:
            if self._type == self._calculated_type and self._end == self._calculated_end and self._len == self._calculated_len:
                self._confidence = 1

            self._sequence = self._alt if self._type == 'INS' or self._calculated_type == 'INS' else (
                self._ref if self._type == 'DEL' or self._calculated_type == 'DEL' else '')

    @property
    def chrom(self):
        return self._chrom

    @property
    def pos(self):
        return self._pos

    @property
    def end(self):
        if not self.indel:
            if self._end:
                return_value = self._end
            else:
                if hasattr(self, '_calculated_end'):
                    return_value = self._pos + self._calculated_end
                else:
                    if hasattr(self, '_calculated_len'):
                        return_value = self._pos + self._calculated_len
                    else:
                        return_value = self._pos + 1
        else:
            return_value = self._end if self._end else self._calculated_end
            if self._type == 'DEL' or self._calculated_type == 'DEL':
                return_value = self._calculated_end
            elif self._type == 'INS' or self._calculated_type == 'INS':
                return_value = self._pos + 1 if not self._end else self._end

        return return_value

    @property
    def sv_type(self):
        return self._type if self._type else self._calculated_type

    def __str__(self):
        return f"{self.chrom}:{self.pos}-{self.end}"

    @property
    def indel(self):
        return (not self._type.split(':')[0].replace('<', '') in ["BND", "DUP", "INV"]) and (not "SNP" in self._info)


def load_variants(vcf_file_path):
    variants = dict()
    with open(vcf_file_path, 'r') as vcf_file:
        for line in vcf_file:
            if not line.startswith("#"):
                values = line.strip().split('\t')
                if len(values) != 10:
                    continue

                # TODO: Remove later
                if values[0] not in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
                                     "16", "17", "18", "19", "X"]:
                    continue

                variant = Variant(values[0], int(values[1]), values[2], values[3], values[4], values[5], values[6],
                                  values[7], values[8], values[9], source=vcf_file_path)
                if variant is None:
                    print(variant)

                lst = variants.get(variant.chrom, dict())
                idx = hash(line)
                lst[idx] = variant
                variants[variant.chrom] = lst

    return variants


def write_variants(merged_variants, output_path):
    with open(output_path, "w") as output_file:
        for variant in merged_variants:
            output_file.write(
                f'{variant.chrom}\t{variant.pos}\t.\t{variant._ref}\t{variant._alt}\t.\tPASS\t{variant._info}\t{variant._format}\t{variant._sample}\n')


class Intersection:
    def __init__(self, line):
        values = line.split('\t')

        self._a = Coordinate(values[0], values[1], values[2])
        self._b = Coordinate(values[4], values[5], values[6])

        self._a_hash = int(values[3])
        self._b_hash = int(values[7])

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @property
    def chrom(self):
        return self._a.chrom

    @property
    def pos(self):
        return self.a.pos if self.a.pos < self.b.pos else self.b.pos

    @property
    def end(self):
        return self.a.end if self.a.end > self.b.end else self.b.end

    @property
    def hash1(self):
        return self._a_hash

    @property
    def hash2(self):
        return self._b_hash

    def __str__(self):
        return f"{self.a.chrom}:{self.a.pos}-{self.a.end} <-> {self.b.chrom}:{self.b.pos}-{self.b.end} => {self.chrom}:{self.pos}-{self.end}"


def merge_variants(variants_a, variants_b, window=0):
    merged_variants = []
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix='.bed') as temp_file_a, \
            tempfile.NamedTemporaryFile(mode="w", delete=False, suffix='.bed') as temp_file_b:

        # print(f"Writing temporary file a: {temp_file_a.name}")
        for idx in variants_a:
            variant = variants_a[idx]
            temp_file_a.write(f"{variant.chrom}\t{int(variant.pos) - window}\t{int(variant.end) + window}\t{idx}\n")

        # print(f"Writing temporary file b: {temp_file_b.name}")
        for idx in variants_b:
            variant = variants_b[idx]
            temp_file_b.write(f"{variant.chrom}\t{int(variant.pos) - window}\t{int(variant.end) + window}\t{idx}\n")

    bed_a = BedTool(temp_file_a.name)
    bed_b = BedTool(temp_file_b.name)

    remove_from_b = []

    for line in bed_a.intersect(bed_b, f=0.8, r=True, wo=True):
        intersection = Intersection(str(line))
        variant_a = variants_a.get(intersection.hash1)
        variant_b = variants_b.get(intersection.hash2)

        if variant_a.sv_type == variant_b.sv_type:
            if variant_a.sv_type in ['INS', 'INDEL']:
                if len(variant_a._sequence) > len(variant_b._sequence):
                    percent = len(variant_b._sequence) / len(variant_a._sequence)
                else:
                    percent = len(variant_a._sequence) / len(variant_b._sequence)

                if percent < 0.8:
                    continue

            variant_a._info += f';GASM={variant_a};PILEUP={variant_b}'
            variant_a._pos = variant_a.pos if variant_a.pos < variant_b.pos else variant_b.pos
            variant_a._end = variant_a.end if variant_a.end > variant_b.end else variant_b.end

            remove_from_b.append(intersection.hash2)

    for entry in list(set(remove_from_b)):
        del variants_b[entry]

    for idx in variants_a:
        variant = variants_a[idx]
        merged_variants.append(variant)

    for idx in variants_b:
        variant = variants_b[idx]
        merged_variants.append(variant)

    return merged_variants


class Merger:

    def __init__(self, gasm_file, pileup_file, output_file):
        self._gasm_file = gasm_file
        self._pileup_file = pileup_file
        self._out_file = output_file

        self._gasm_variants = []
        self._pileup_variants = []

    def merge(self):
        merged_variants = []
        logger.info(f"Loading entries from file {self._gasm_file}")
        self._gasm_variants = load_variants(self._gasm_file)
        logger.info(f"Loading entries from file {self._pileup_file}")
        self._pileup_variants = load_variants(self._pileup_file)

        for chrom in self._gasm_variants.keys():
            logger.info(f"Merging chromosome {chrom}")
            merged_variants.extend(merge_variants(self._gasm_variants.get(chrom), self._pileup_variants.get(chrom), 5))

        logger.info(f"Writing result to {self._out_file}")
        write_variants(merged_variants, self._out_file)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge variants from mpileup and minigraph + svim-asm",
        epilog="Example: python merge.py --gasm_file=svim.vcf --pileup_file=pileup_del.vcf",
    )
    parser.add_argument(
        "--gasm_file",
        type=Path,
        help="Input VCF (minigraph + svim-asm)",
    )
    parser.add_argument(
        "--pileup_file",
        type=Path,
        help="Input VCF (pileup)",
    )
    parser.add_argument(
        "--out",
        type=Path,
        help="Output VCF",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.gasm_file.is_file():
        logger.error(f"The given input file {args.gasm_file} was not found!")
        sys.exit(2)

    if not args.pileup_file.is_file():
        logger.error(f"The given input file {args.pileup_file} was not found!")
        sys.exit(2)

    Merger(args.gasm_file, args.pileup_file, args.out).merge()


if __name__ == "__main__":
    sys.exit(main())
