import os
import unittest

from samsum import testing_utils as utils


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        from samsum import commands
        self.test_fa = utils.get_test_data("fasta_test.fa")
        self.test_sam = utils.get_test_data("samsum_test_2.sam")
        self.test_ref_fa = utils.get_test_data("samsum_test_2.fasta")
        self.output_table = os.path.join(utils.get_project_root(), "tests/samsum_table.csv")

        self.ref_seq_abundances = commands.ref_sequence_abundances(aln_file=self.test_sam, seq_file=self.test_ref_fa,
                                                                   min_aln=10, p_cov=0, map_qual=0)
        return

    def test_write_summary_table(self):
        """ Ensure the header and number of columns in the output table is correct """
        from samsum import file_parsers as ss_fp
        curr_table_header = ["QueryName", "RefSequence", "ProportionCovered", "Coverage", "Fragments", "FPKM", "TPM"]
        ss_fp.write_summary_table(self.ref_seq_abundances, self.output_table, "pytest", 0)
        with open(self.output_table) as table_handler:
            header_fields = table_handler.readline().strip().split(',')
            data_lines = []
            line = table_handler.readline()
            while line:
                data_lines.append(line.strip())
                line = table_handler.readline()

        self.assertEqual(header_fields, curr_table_header)
        for line in data_lines:
            self.assertEqual(len(line.split(',')), len(curr_table_header))
        return

    def test_ref_sequence_length(self) -> None:
        """
        Ensure the RefSequence.rightmost doesn't exceed its length from alignment_dat_example.load_sam()

        :return:
        """
        from samsum import file_parsers as ss_fp
        ref_seq_lengths = ss_fp.fasta_seq_lengths(fasta_file=self.test_ref_fa)
        self.assertEqual(277, len(ref_seq_lengths))

        # Parse the alignments and return the strings of reads mapped to each reference sequence
        mapped_dict = ss_fp.sam_parser_ext(self.test_sam, True)
        unmapped_aln = mapped_dict.pop("UNMAPPED").pop()
        self.assertEqual(4890.5, unmapped_aln.weight)

        # Convert the alignment strings returned by the sam_parser_ext into ss_class.AlignmentDat instances
        for matches in mapped_dict.values():  # type: list
            for match in matches:
                self.assertTrue(match.end < ref_seq_lengths[match.subject])
        return

    def test_fasta_reader(self) -> None:
        from samsum import file_parsers as ss_fp
        ref_seq_lengths = ss_fp.fasta_seq_lengths(fasta_file=self.test_fa)
        self.assertEqual({"contig 1": 16, "contig 2": 25}, ref_seq_lengths)
        return


if __name__ == '__main__':
    unittest.main()
