from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class FastQFilter:
    def __init__(
        self, input_file, output_file, min_length=0, min_quality=0, min_gc=0, max_gc=100
    ):
        self.input_file = input_file
        self.output_file = output_file
        self.min_length = min_length
        self.min_quality = min_quality
        self.min_gc = min_gc
        self.max_gc = max_gc

    def filter_fastq(self):
        with open(self.output_file, "w") as output_handle:
            for record in SeqIO.parse(self.input_file, "fastq"):
                if self._passes_filter(record):
                    SeqIO.write(record, output_handle, "fastq")

    def _passes_filter(self, record):
        if len(record.seq) < self.min_length:
            return False

        if min(record.letter_annotations["phred_quality"]) < self.min_quality:
            return False

        gc_content = gc_fraction(record.seq)
        if gc_content < self.min_gc or gc_content > self.max_gc:
            return False

        return True
