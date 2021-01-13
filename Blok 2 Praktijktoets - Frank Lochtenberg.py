# Made by Frank Lochtenberg

import re
import matplotlib.pyplot as plt
import tkinter
from tkinter import messagebox


def read_fasta_header():
    """ Reads the fasta file and puts the amino acids headers in a list.

    :return: list of headers of amino acids.
    """
    fasta = open("Caenorhabditis_elegans.cds_pep.all.fa", "r")
    headers = []
    for line in fasta:
        if ">" in line:
            headers.append(line)
    return headers


def read_fasta_seq():
    """ Reads the fasta-file and puts the amino acid sequences in a list

    :return: list of amino acids sequences.
    """
    fasta = open("Caenorhabditis_elegans.cds_pep.all.fa", "r")
    t_seq = []
    d_seq = []
    for line in fasta:
        line = line.strip()
        if ">" in line:
            if t_seq != []:
                d_seq.append("".join(t_seq))
                t_seq = []
        else:
            t_seq.append(line)
    d_seq.append("".join(t_seq))
    return d_seq


def counting_exons():
    """Counts the frequency of exons in the gff3 file.

    :return: An integer which is the frequency of exons in the gff3 file
    """
    gff3 = open("Caenorhabditis_elegans.gff3", "r")
    exon_count = 0
    for line in gff3:
        if "x" in line[11:14]:
            exon_count += 1
    return exon_count


def counting_cds():
    """Counts the frequency CDS in the gff3 file.

    :return: An integer which is the frequency of CDS in the gff3 file.
    """
    gff3 = open("Caenorhabditis_elegans.gff3", "r")
    cds_count = 0
    for line in gff3:
        if "D" in line[11:14]:
            cds_count += 1
    return cds_count


def counting_mrna():
    """Counts the frequency of mRNA in the gff3 file.

    :return: An integer which is the frequency of mRNA in the gff3 file
    """
    gff3 = open("Caenorhabditis_elegans.gff3", "r")
    mrna_count = 0
    for line in gff3:
        if "m" in line[11:14]:
            mrna_count += 1
    return mrna_count


def display_count(exon_count, cds_count, mrna_count):
    """Displays the frequency of exons, CDS, mRNA in the gff3 file.

    :param exon_count: An integer which is the frequency of exons in
    the gff3 file.
    :param cds_count: An integer which is the frequency of CDS in
    the gff3 file.
    :param mrna_count: An integer which is the frequency of mRNA in
    the gff3 file.
    :return: Nothing
    """
    print(exon_count, ",", cds_count, "&", mrna_count)


class concencus_check:
    """Checks if the is a concencus pattern in the amino acids sequences

    """

    def set_headers(self, headers):
        self.__headers = headers

    def get_headers(self):
        return self.__headers

    def set_aa_seq(self, seq):
        self.__aa_seq = seq

    def get_aa_seq(self):
        return self.__aa_seq

    def set_concensus(self):
        for aa_seq in self.get_aa_seq():
            ccc = re.search(r"C[A-Z]{2}C[A-Z]{2}C[A-Z]{5}C[A-Z]{2}"
                            r"C[A-Z]{2}", aa_seq)
            if ccc:
                self.__concencus = ccc
            else:
                print("No!")

    def get_concencus(self):
        return self.__concencus

    def print_concencus(self):
        print(self.get_concencus())

class GUI:
    # Maakt een main window.
    self.main_window = tkinter.Tk()
    self.main_window.title("km to miles converter")


def main():
    headers = read_fasta_header()
    d_seq = read_fasta_seq()
    cccc = concencus_check()
    cccc.set_headers(headers)
    cccc.set_aa_seq(d_seq)
    exon_count = counting_exons()
    cds_count = counting_cds()
    mrna_count = counting_mrna()
    display_count(exon_count, cds_count, mrna_count)


main()
