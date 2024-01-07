"""
Author  : D. R. R. Wijewardene      S14245
Date    : 12/01/2022
Project2 - Sequence Aligner
Input   : A FASTA file containing multiple DNA, RNA, or amino acid sequences
Output  :
    Type of FASTA sequence
    Global pairwise similarity score
    Scatter plot of pairwise similarity
    BLAST search results of sequences in FASTA file
    Sequence obtained from GenBank
    Melting temperature of nucleotide sequence
    Isoelectric point of amino acid sequence
    Aromaticity of amino acid sequence
"""

import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO, pairwise2, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils.ProtParam import ProteinAnalysis as pa


class SeqFile:

    # A constructor method to create sequence objects
    def __init__(self, seq_file):
        self.records = list(SeqIO.parse(seq_file, "fasta"))

    # 1 - method to output the type
    @staticmethod
    def get_seq_type(sequence):
        seq = str(sequence)
        # default
        sq_type = "Amino Acid"
        # Declare DNA and mRNA base patterns
        dna_bases = re.compile("^[ATGC]*$")
        mrna_bases = re.compile("^[AUGC]*$")
        # Assign the sequence type
        if dna_bases.match(seq):
            sq_type = "DNA"
        elif mrna_bases.match(seq):
            sq_type = "RNA"
        # Return the sequence type
        return sq_type

    # 2 - method to calculate the global pairwise similarity score
    def get_gps_score(self):
        seq1_list = []
        seq2_list = []
        score_list = []
        # obtain sequence pairs to calculate the score
        for seq1 in self.records:
            seq1_type = self.get_seq_type(seq1.seq)
            for seq2 in self.records:
                seq2_type = self.get_seq_type(seq2.seq)
                if not seq1.seq == seq2.seq and seq1_type == seq2_type and seq2.id not in seq1_list:
                    # calculate score
                    alignment_score = pairwise2.align.globalxx(seq1.seq, seq2.seq, score_only=True)
                    seq1_list.append(seq1.id)
                    seq2_list.append(seq2.id)
                    score_list.append(alignment_score)
        # create dataframe with sequence pairs and their alignment score
        gpss = pd.DataFrame({'Sequence 1': seq1_list, 'Sequence 2': seq2_list, 'Score': score_list})
        # write the output in a text file
        # noinspection PyTypeChecker
        np.savetxt('Outputs/global_pss.txt', gpss.values, fmt='%s', delimiter='\t',
                   header='Sequence 1\tSequence 2\tScore')
        return gpss

    # 3 - method to generate a scatter plot to represent pairwise similarity
    def scatter_plot(self, seq1, seq2):
        seq_1 = seq1.seq
        seq_2 = seq2.seq
        fig, axes = plt.subplots(figsize=(8, 6))
        seq1_type = self.get_seq_type(seq_1)
        seq2_type = self.get_seq_type(seq_2)
        if seq1_type == seq2_type:
            window = 7
            seq_one = seq_1.upper()
            seq_two = seq_2.upper()
            # create dictionaries mapping the window-sized sub-sequences to locations
            dict_one = {}
            dict_two = {}
            for (seq, section_dict) in [
                (seq_one, dict_one),
                (seq_two, dict_two)
            ]:
                for i in range(len(seq) - window):
                    section = seq[i: i + window]
                    try:
                        section_dict[section].append(i)
                    except KeyError:
                        section_dict[section] = [i]
            # find any sub-sequences found in both sequences
            matches = set(dict_one).intersection(dict_two)
            # create lists of x and y co-ordinates for scatter plot
            x = []
            y = []
            for section in matches:
                for i in dict_one[section]:
                    for j in dict_two[section]:
                        x.append(i)
                        y.append(j)
            # create dataframe with the data for the plot
            x_name = seq1.id
            y_name = seq2.id
            df = pd.DataFrame({x_name: x, y_name: y})
            # create plot
            if not df.empty:
                sns.scatterplot(data=df, x=x_name, y=y_name, ax=axes, color="purple")
                axes.set_xlim(0, len(seq1) - window)
                axes.set_ylim(0, len(seq2) - window)
                axes.set_title("Scatter Plot showing Pairwise Similarity between %s and %s\n(allowing no mis-matches)"
                               % (x_name, y_name))
                plt.savefig("Outputs/scatter_plot_%s_and_%s.png" % (x_name, y_name))
                plt.show()
        else:
            print("Pairwise Similarity Scatter Plot - Not Applicable for %s and %s" % (seq1.id, seq2.id))
        return

    # 4 - method to run BLAST searches
    def seq_blast(self):
        for i in range(len(self.records)):
            seq = str(self.records[i].seq)
            seq_type = self.get_seq_type(seq)
            if seq_type == "Amino Acid":
                # Run the web-based nucleotide BLAST program
                result_handle = NCBIWWW.qblast("blastp", "nr", seq)
                # Parse through the BLAST hits in the file
                blast_records = NCBIXML.parse(result_handle)
            else:
                # Run the web-based nucleotide BLAST program
                result_handle = NCBIWWW.qblast("blastn", "nt", seq)
                # Parse through the BLAST hits in the file
                blast_records = NCBIXML.parse(result_handle)
            # Save blast results in separate document
            with open("Outputs/blast_results.xml", "w") as out_handle:
                out_handle.write(result_handle.read())
            result_handle.close()
            return blast_records

    # 5 - method to find a sequence from GenBank when its accession number is given
    @staticmethod
    def get_from_genbank(acc_no):
        Entrez.email = "roshw1998@gmail.com"
        # get the sequence
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=acc_no)
        # save in a fasta file
        seq_record = SeqIO.read(handle, "genbank")
        new_record_file = open("Outputs/%s_record.fasta" % seq_record.id, "w")
        new_record_file.write(">%s\n%s" % (seq_record.description, seq_record.seq))
        handle.close()
        return

    # 6 - method to calculate the melting temperature of a given nucleotide sequence
    def melting_temp(self, seq):
        seq_type = self.get_seq_type(seq)
        if seq_type == "DNA" or seq_type == "RNA":
            mel_temp = mt.Tm_NN(seq)
            print("Melting Temperature : %0.2f \N{DEGREE SIGN}C" % mel_temp)
            return mel_temp
        else:
            print("Melting Temperature : Not applicable for %s" % seq.id)
            return

    # 7 - method to calculate the isoelectric point of a given amino acid sequence
    def isoelectric(self, sequence):
        seq = str(sequence)
        seq_type = self.get_seq_type(seq)
        if seq_type == "Amino Acid":
            ie_point = pa(seq).isoelectric_point()
            print("Isoelectric Point : %0.2f" % ie_point)
            return ie_point
        else:
            print("Isoelectric Point : Not applicable for %s" % sequence.id)
            return

    # 8 - method to calculate the aromaticity of a given amino acid sequence
    def aromaticity(self, sequence):
        seq = str(sequence)
        seq_type = self.get_seq_type(seq)
        if seq_type == "Amino Acid":
            aroma = pa(seq).aromaticity()
            print("Aromaticity : %0.2f" % aroma)
            return aroma
        else:
            print("Aromaticity : Not applicable for %s" % sequence.id)
            return
