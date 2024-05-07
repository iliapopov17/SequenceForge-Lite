import unittest
from sequence_forge import fastq_filter, DNASequence, RNASequence, AminoAcidSequence
from bio_files_processor import (
    convert_multiline_fasta_to_oneline,
    change_fasta_start_pos,
)
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os


class TestSequenceForge(unittest.TestCase):
    """
    Test suite for sequence_forge module, focusing on functionality, error handling, and file operations.
    """

    def test_fastq_filter_correctness(self):
        """
        Test if fastq_filter correctly filters records based on specified criteria.
        """
        # Setup: Create a dummy FASTQ file with mock data
        records = [
            SeqRecord(
                Seq("AGCT" * 100),
                id="high_gc",
                description="high_gc_content",
                letter_annotations={"phred_quality": [40] * 400},
            ),
            SeqRecord(
                Seq("ATAT" * 100),
                id="low_gc",
                description="low_gc_content",
                letter_annotations={"phred_quality": [20] * 400},
            ),
        ]
        input_path = "test_input.fastq"
        output_path = "test_output.fastq"
        SeqIO.write(records, input_path, "fastq")

        # Action: Filter the records
        fastq_filter(
            input_fastq=input_path,
            output_fastq=output_path,
            gc_bound=(50, 100),
            length_bound=(300, 500),
            quality_threshold=30,
        )

        # Check: Verify output file has only the high GC content record
        output_records = list(SeqIO.parse(output_path, "fastq"))
        self.assertEqual(len(output_records), 1)
        self.assertEqual(output_records[0].id, "high_gc")

        # Cleanup
        os.remove(input_path)
        os.remove(output_path)

    def test_fastq_filter_error_handling(self):
        """
        Test if fastq_filter handles non-existent input files gracefully.
        """
        with self.assertRaises(FileNotFoundError):
            fastq_filter(
                input_fastq="non_existent.fastq", output_fastq="irrelevant.fastq"
            )

    def test_dna_sequence_transcription(self):
        """
        Test if DNASequence transcribes correctly to RNASequence.
        """
        dna = DNASequence("ATGC")
        rna = dna.transcribe()
        self.assertEqual(str(rna), "UACG")

    def test_amino_acid_sequence_validation(self):
        """
        Test if AminoAcidSequence checks alphabet validity correctly.
        """
        aa_seq = AminoAcidSequence("ACDEFGHIKLMNPQRSTVWY")
        self.assertTrue(aa_seq.check_alphabet())

    def test_amino_acid_sequence_invalid_input(self):
        """
        Test error handling for invalid amino acid sequences.
        """
        aa_seq = AminoAcidSequence("ACDEFGHIKLMNPQRSTVWYZ")
        self.assertFalse(aa_seq.check_alphabet())

    def test_file_operations_amino_acid_sequence(self):
        """
        Test reading and writing amino acid sequences to a file.
        """
        aa_seq = AminoAcidSequence("ACDEFGHIKLMNPQRSTVWY")
        with open("temp_aa_seq.txt", "w") as file:
            file.write(str(aa_seq))
        with open("temp_aa_seq.txt", "r") as file:
            content = file.read()
        self.assertEqual(content, "ACDEFGHIKLMNPQRSTVWY")
        os.remove("temp_aa_seq.txt")


class TestBioFilesProcessor(unittest.TestCase):
    """
    Test suite for bio_files_processor module, specifically testing FASTA file manipulation functions.
    """

    def test_convert_multiline_fasta_to_oneline(self):
        """
        Test that a multiline FASTA file is correctly converted to a single-line format.
        """
        # Setup: Create a multiline FASTA file
        multiline_fasta = "test_multiline.fasta"
        oneline_fasta = "test_oneline.fasta"
        with open(multiline_fasta, "w") as f:
            f.write(">seq1\nATGC\nATGC\n>seq2\nGCTA\nGCTA\n")

        # Action: Convert it to a single-line FASTA
        convert_multiline_fasta_to_oneline(multiline_fasta, oneline_fasta)

        # Check: Verify the contents of the output file
        with open(oneline_fasta, "r") as f:
            lines = f.readlines()
        expected_lines = [">seq1\n", "ATGCATGC\n", ">seq2\n", "GCTAGCTA\n"]
        self.assertEqual(lines, expected_lines)

        # Cleanup
        os.remove(multiline_fasta)
        os.remove(oneline_fasta)

    def test_change_fasta_start_pos(self):
        """
        Test shifting the start position of sequences in a FASTA file.
        """
        # Setup: Create a single-line FASTA file
        fasta_file = "test_fasta.fasta"
        shifted_fasta = "test_shifted.fasta"
        with open(fasta_file, "w") as f:
            f.write(">seq1\nATGCATGC\n")

        # Action: Shift the sequence start position
        change_fasta_start_pos(fasta_file, 4, shifted_fasta)

        # Check: Verify the contents of the shifted file
        with open(shifted_fasta, "r") as f:
            lines = f.readlines()
        expected_lines = [">seq1\n", "ATGCATGC\n"]  # Expected shifted sequence
        self.assertEqual(lines, expected_lines)

        # Cleanup
        os.remove(fasta_file)
        os.remove(shifted_fasta)


if __name__ == "__main__":
    unittest.main()
