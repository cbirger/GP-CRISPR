# -*- coding: utf-8 -*-
"""
Trims FASTQ-formatted reads from a CRISPR pooled-screen experiment so only sgRNA sequences remain.

The fuzzy search (tolerate up to two mismatches) for prefix using regular expression search
was derived from Galaxy tool seq_primer_clip, available via the Galaxy tool Shed at
http://toolshed.g2.bx.psu.edu/view/peterjc/seq_primer_clip. 

"""

import argparse
import sys
import os
import gzip
import re

from Bio import SeqIO, Seq

#---------------------
#   Classes
#---------------------

class SgRNAExtractor(object):

    def __init__(self, reads_file, reverse_reads, prefix_sequence,
                 prefix_mismatches, sgRNA_length, max_reads, output_basename):
        """
        Args:
          reads_file: Filename or input file handle of fastq or compressed (gz)
            fastq file.
          reverse_reads: Boolean variable set to true of reads_file is reverse
            reads file from dual CRISPR experiment.
          prefix_sequence: String containing nucleotide sequence that marks
            start of the sgRNA.
          prefix_mismatches: Maximum number of mismatches allowed when
            matching prefix to read sequence.
          sgRNA_length: Length of sgRNA
          max_reads: The maximum number of read records to process.
          output_basename: The basename to use for two output files.
          
        """
        self.reads_file = reads_file
        self.reverse_reads = reverse_reads
        self.prefix_sequence = prefix_sequence
        self.max_mismatches = prefix_mismatches
        self.sgRNA_length = sgRNA_length
        self.max_reads = max_reads
        self.output_basename = output_basename

        self.record_count = 0
        self.clipped_count = 0
        self.unclipped_count = 0

        filename, fileExtension = os.path.splitext(self.reads_file)
        if fileExtension == ".gz":
            self.gzipped = True
        else:
            self.gzipped = False

        if self.gzipped:
            filename, fileExtension = os.path.splitext(self.output_basename)
            if fileExtension == "fastq":
                self.output_basename = filename

        ambiguous_dna_values = {
            "A": "A",
            "C": "C",
            "G": "G",
            "T": "T",
            "M": "ACM",
            "R": "AGR",
            "W": "ATW",
            "S": "CGS",
            "Y": "CTY",
            "K": "GTK",
            "V": "ACGMRSV",
            "H": "ACTMWYH",
            "D": "AGTRWKD",
            "B": "CGTSYKB",
            "X": ".", #faster than [GATCMRWSYKVVHDBXN] or even [GATC]
            "N": ".",
            }

        self.ambiguous_dna_re = {}
        for letter, values in ambiguous_dna_values.iteritems():
            if len(values) == 1:
                self.ambiguous_dna_re[letter] = values
            else:
                self.ambiguous_dna_re[letter] = "[%s]" % values

    def extract_sgRNAs(self):
        """
        Extracts sgRNA sequences from FASTQ records and writes them to a new FASTQ file.
        """
        
        if not self.gzipped:
            original_records = SeqIO.parse(self.reads_file, "fastq")
            trimmed_records = self.__clip_reads(original_records)
            count = SeqIO.write(trimmed_records, self.output_basename + ".fastq", "fastq")
            return count
        else:
            handle = gzip.open(self.reads_file)
            original_records = SeqIO.parse(handle, "fastq")
            trimmed_records = self.__clip_reads(original_records)
            count = SeqIO.write(trimmed_records, self.output_basename + ".fastq", "fastq")
            handle.close()
            return count

    def __clip_reads(self, input_record_iterator):

        re_pattern = self.__load_prefix_as_re()
        self.record_count = 0
        self.clipped_count = 0
        self.unclipped_count = 0

        while self.record_count != self.max_reads:
            record = input_record_iterator.next()
            self.record_count += 1
            # use regular expression search to find prefix match
            result = re_pattern.search(str(record.seq.upper()))
            if result:
                # when searching for prefix in forward read, first sgRNA is directly after prefix
                # when searching for prefix in reverse read, second sgRNA is directly after prefix
                cut = result.end()
                if len(record) - cut >= self.sgRNA_length:
                    self.clipped_count += 1
                    yield record[cut:cut+self.sgRNA_length]
                else:
                    # while found pattern match, the number of bases in sequence directly after
                    # matching prefix less than sgRNA length
                    # in this case, we just yeild first sgRNA_length bases in order to keep pairing
                    # they won't align to ref sgRNA sequences
                    self.unclipped_count += 1
                    yield record[0:self.sgRNA_length]
            else:
                self.unclipped_count += 1
                # if no pattern match, just yield sgRNA_length bases
                # we do this to keep pairing - they won't align to known ref sgRNA sequences
                yield record[0:self.sgRNA_length]

    def __load_prefix_as_re(self):
        re_patterns = set()
        prefix_seq = Seq.Seq(self.prefix_sequence)
        if self.reverse_reads:
            prefix_seq = prefix_seq.reverse_complement()
        for pattern in self.__make_reg_ex_mm(prefix_seq, self.max_mismatches):
            re_patterns.add(pattern)
    
        # Use set to avoid duplicates, sort to have longest first
        # (so more specific primers found before less specific ones)
        re_patterns = sorted(set(re_patterns), key=lambda p: -len(p))
        return re.compile("|".join(re_patterns)) # make one monster re!

    def __make_reg_ex_mm(self, prefix_seq, mm):

        if mm > 2:
            raise NotImplementedError("At most 2 mismatches allowed!")

        seq = prefix_seq.upper()
        yield self.__make_reg_ex(seq)
        
        for i in range(1,mm+1):
            # Missing first/last i bases at very start/end of sequence
            for reg in self.__make_reg_ex_mm(seq[i:],  mm-i):
                yield "^" + reg
            for reg in self.__make_reg_ex_mm(seq[:-i], mm-i):
                yield "$" + reg
                
        if mm >= 1:
            for i,letter in enumerate(seq):
                # duplicate patterns will be removed later
                # if letter not in "NX":
                pattern = seq[:i] + "N" + seq[i+1:]
                assert len(pattern) == len(seq), "Len %s is %i, len %s is %i" \
                   % (pattern, len(pattern), seq, len(seq))
                yield self.__make_reg_ex(pattern)

        if mm >=2:
            for i,letter in enumerate(seq):

                for k,letter in enumerate(seq[i+1:]):
                    pattern = seq[:i] + "N" + seq[i+1:i+1+k] + "N" + seq[i+k+2:]
                    assert len(pattern) == len(seq), "Len %s is %i, len %s is %i" \
                       % (pattern, len(pattern), seq, len(seq))
                    yield self.__make_reg_ex(pattern)

    def __make_reg_ex(self, ambiguous_sequence):
        return "".join(self.ambiguous_dna_re[letter] for letter in ambiguous_sequence)

                                    
    def write_summary_stats(self):
        """
        Writes input parameters and trimming statistics to a text file.
        """

        handle = open(self.output_basename + "_SUMMARY.txt", 'w')
        

        # input parameters
        handle.write("INPUT PARAMETERS\n\n")
        handle.write("reads_file = %s\n" % self.reads_file)
        handle.write("prefix = %s\n" % self.prefix_sequence)
        handle.write("output_basename = %s\n" % self.output_basename)
        handle.write("max_reads = %i\n" % self.max_reads)
        handle.write("max_mismatches = %i\n" % self.max_mismatches)
        handle.write("reverse_reads = %s\n" % self.reverse_reads)

        handle.write("\nOUTPUT STATISTICS\n\n")
        handle.write("read count: %i\n" % self.record_count)
        handle.write("clipped reads: %i\n" % self.clipped_count)
        handle.write("unclipped reads: %i\n" % self.unclipped_count)

        handle.close()


