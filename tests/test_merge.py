from unittest import TestCase
from unittest.mock import Mock

from bin.merge import Variant


class TestVariant(TestCase):

    def test__get_precise_len_and_end(self):
        variant = Mock(spec=Variant,
                       _pos=10,
                       _ref='A',
                       _alt='ACTCAAATCAGGGAGTAGAAAACACTGCCCTTTTTCTCCCTACCTATTCAAT',
                       _info='SVTYPE=INS;END=3094520;SVLEN=99')
        sv_type, sv_len, end = Variant._get_precise_type_len_end(variant)

        self.assertEqual('INS', sv_type)
        self.assertEqual(51, sv_len)
        self.assertEqual(61, end)

    def test__process_info_complete(self):
        variant = Mock(spec=Variant,
                       _info='SVTYPE=INS;END=3094510;SVLEN=51')

        sv_type, sv_len, end, = Variant._get_info_type_len_end(variant)
        self.assertEqual('INS', sv_type)
        self.assertEqual(51, sv_len)
        self.assertEqual(3094510, end)

    def test__process_info_incomplete_no_end(self):
        variant = Mock(spec=Variant,
                       _info='SVTYPE=INS;SVLEN=51')

        sv_type, sv_len, end, = Variant._get_info_type_len_end(variant)
        self.assertEqual('INS', sv_type)
        self.assertEqual(51, sv_len)
        self.assertEqual(0, end)

    def test__process_info_incomplete_no_type(self):
        variant = Mock(spec=Variant,
                       _info='SVLEN=51')

        sv_type, sv_len, end, = Variant._get_info_type_len_end(variant)
        self.assertEqual('', sv_type)
        self.assertEqual(51, sv_len)
        self.assertEqual(0, end)

    def test__get_target_sequence_ins(self):
        variant = Mock(spec=Variant,
                       _ref='G',
                       _alt='GGATA',
                       _type='INS')
        target_sequence = Variant._get_target_sequence(variant)
        self.assertEqual('GGATA', target_sequence)

    def test__get_target_sequence_del(self):
        variant = Mock(spec=Variant,
                       _ref='GATA',
                       _alt='A',
                       _type='DEL')
        target_sequence = Variant._get_target_sequence(variant)
        self.assertEqual('GATA', target_sequence)