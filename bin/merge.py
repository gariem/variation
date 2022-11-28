import argparse
import logging
import sys
from pathlib import Path

logger = logging.getLogger()


class Variant:
    def __init__(self, chrom, pos, identifier, ref, alt, quality, info, format_value, sample, sv_type='', caller='',
                 precise=True):
        self._chrom = chrom
        self._pos = pos
        self._id = identifier
        self._ref = ref
        self._alt = alt
        self._quality = quality
        self._info = info
        self._format = format_value
        self._sample = sample
        self._caller = caller
        self._precise = precise

        self._type, self._len, self._end = self._get_precise_type_len_end() if self._precise \
            else self._get_info_type_len_end()

        self._sequence = self._get_target_sequence()

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, Variant):
            return self._chrom == other._chrom and self._pos == other._pos
        return False

    def _get_precise_type_len_end(self):
        sv_len = 0
        end = 0

        if self._ref and self._alt and self._pos:
            sv_len = len(self._alt) - len(self._ref)
            sv_type = 'INS' if sv_len > 0 else ('DEL' if sv_len < 0 else '')
            end = self._pos + abs(sv_len)

        return sv_type, sv_len, end

    def _get_info_type_len_end(self):
        sv_type = ''
        sv_len = 0
        end = 0

        info_list = self._info.split(';')

        for info in info_list:
            if 'SVTYPE=' in info:
                sv_type = info.replace('SVTYPE=', '')
            if 'SVLEN=' in info:
                sv_len = int(info.replace('SVLEN=', ''))
                continue
            if 'END=' in info:
                end = int(info.replace('END=', ''))
                continue

        return sv_type, sv_len, end

    def _get_target_sequence(self):
        return self._ref if self._type == 'DEL' else self._alt


class Merger:
    def __init__(self, gasm_file, mpileup_file):
        self._gasm_file = gasm_file
        self._mpileup_file = mpileup_file

    def merge(self):
        print("Hello world!!")


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge variants from mpileup and minigraph + svim-asm",
        epilog="Example: python merge.py --gasm_file=svim.vcf --pileup_file=pileup.vcf",
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

    Merger(args.gasm_file, args.pileup_file).merge()


if __name__ == "__main__":
    sys.exit(main())
