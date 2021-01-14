# Made by Frank Lochtenberg
# Version 4

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


def counting_total():
    """Counts the total lines in the gff3 file.

    :return: An integer which is the total of lines in the gff3 file.
    """
    gff3 = open("Caenorhabditis_elegans.gff3", "r")
    total_count = 0
    for line in gff3:
        if line != "":
            total_count += 1
    return total_count


def calculating_other(exon_count, cds_count, mrna_count, total_count):
    """ Calculates the frequency of other items than exons, CDS & mRNA
    in the gff3 file.

    :param exon_count: An integer which is the frequency of exons in
    the gff3 file.
    :param cds_count: An integer which is the frequency of CDS in
    the gff3 file.
    :param mrna_count: An integer which is the frequency of mRNA in
    the gff3 file.
    :param total_count: An integer which is the total of lines in
    the gff3 file.
    :return: An integer which is the frequency of other items than exons
    , CDS & mRNA in the gff3 file.
    """
    other_count = total_count - exon_count - cds_count - mrna_count
    return other_count


class ConcencusCheck:
    """Checks if the is a concencus pattern in the amino acids sequences

    """

    def __init__(self, headers, seq):
        self.set_headers(headers)
        self.set_aa_seq(seq)
        self.set_concensus()
        self.print_concencus()

    def set_headers(self, headers):
        self.headers = headers

    def get_headers(self):
        return self.headers

    def set_aa_seq(self, seq):
        self.aa_seq = seq

    def get_aa_seq(self):
        return self.aa_seq

    def set_concensus(self):
        nccc = 0
        cccl = []
        for aa_seq in self.get_aa_seq():
            ccc = re.search(r"C[A-Z]{2}C[A-Z]{2}C[A-Z]{5}C[A-Z]{2}"
                            r"C[A-Z]{2}C", aa_seq)
            if ccc:
                cccl.append(ccc.group())
                self.concencus = cccl
            else:
                nccc += 1
        return nccc

    def get_concencus(self):
        return self.concencus

    def print_concencus(self):
        lccc = self.get_concencus()
        print(self.get_concencus())
        return lccc

class GUI:
    def __init__(self, exon_count, cds_count, mrna_count, total_count,
                 other_count):
        self.exons_count = exon_count
        self.cds_count = cds_count
        self.mrna_count = mrna_count
        self.total_count = total_count
        self.other_count = other_count
        headers = read_fasta_header()
        seq = read_fasta_seq()
        zfc = ConcencusCheck(headers, seq).get_concencus()
        self.zfc = zfc
        self.show_plot = plt.show()

        # Makes the main window
        self.main_window = tkinter.Tk()
        self.main_window.title("Information over the sequences")

        # Makes frames
        self.intro_frame = tkinter.Frame(self.main_window)
        self.info_buttons_frame = tkinter.Frame(self.main_window)
        self.graph_buttons_frame = tkinter.Frame(self.main_window)
        self.quit_frame = tkinter.Frame(self.main_window)

        # Makes a label for the introduction sentence
        self.intro_label = tkinter.Label(self.intro_frame,
                                         text="Click the buttons for "
                                              "information!")

        # Places the label for the introduction sentence
        self.intro_label.pack()

        # Makes the buttons for info
        self.exon_button = tkinter.Button(self.info_buttons_frame,
                                          text="Exons Info",
                                          command=self.exons)
        self.cds_button = tkinter.Button(self.info_buttons_frame,
                                         text="CDS Info",
                                         command=self.cds)
        self.mrna_button = tkinter.Button(self.info_buttons_frame,
                                          text="mRNA Info",
                                          command=self.mrna)
        self.other_button = tkinter.Button(self.info_buttons_frame,
                                           text="Other items Info",
                                           command=self.other)
        self.total_button = tkinter.Button(self.info_buttons_frame,
                                           text="Total of items Info",
                                           command=self.total)
        self.concencus_button = tkinter.Button(self.info_buttons_frame,
                                               text="Zinc Finger "
                                                    "Concencus Info",
                                               command=self.concencus)

        # Places the buttons for info
        self.exon_button.pack(side="left")
        self.cds_button.pack(side="left")
        self.mrna_button.pack(side="left")
        self.other_button.pack(side="left")
        self.total_button.pack(side="left")
        self.concencus_button.pack(side="left")

        # Makes the graph buttons
        self.pie_button = tkinter.Button(self.graph_buttons_frame,
                                         text="Show Pie Diagram",
                                         command=self.pie_diagram)

        # Places the graph buttons
        self.pie_button.pack(side="left")

        # Makes the quit button
        self.quit_button = tkinter.Button(self.quit_frame,
                                          text="Quit",
                                          command=self.main_window.
                                          destroy)
        # Places the quit button
        self.quit_button.pack(side="right")

        # Places the frames in the GUI.
        self.intro_frame.pack()
        self.info_buttons_frame.pack()
        self.graph_buttons_frame.pack()
        self.quit_frame.pack()

        # Makes the main window visible.
        tkinter.mainloop()

    def exons(self):
        tkinter.messagebox.showinfo("Exons Info",
                                    "There are " +
                                    str(self.exons_count) +
                                    " exon sequences in the "
                                    "C. Elegans.")

    def cds(self):
        tkinter.messagebox.showinfo("CDS Info",
                                    "There are " +
                                    str(self.cds_count) +
                                    " CDS sequences in the C. Elegans.")

    def mrna(self):
        tkinter.messagebox.showinfo("mRNA Info",
                                    "There are " +
                                    str(self.mrna_count) +
                                    " mRNA sequences in the "
                                    "C. Elegans.")

    def other(self):
        tkinter.messagebox.showinfo("Other items Info",
                                    "There are " +
                                    str(self.other_count) +
                                    " of other items than exon, CDS & "
                                    "mRNA sequences in the C. Elegans.")

    def total(self):
        tkinter.messagebox.showinfo("Total of items Info",
                                    "There are " +
                                    str(self.total_count) +
                                    " of total sequences in the "
                                    "C. Elegans.")

    def concencus(self):
        zinc_finger_count = 0
        for _ in self.zfc:
            zinc_finger_count += 1
        tkinter.messagebox.showinfo("Zinc Finger Concencus Info",
                                    "There are " +
                                    str(zinc_finger_count) +
                                    " of the Zinc Finger Concencus "
                                    "constated in the C. Elegans.")

    def pie_diagram(self):
        values = [self.exons_count, self.cds_count, self.mrna_count,
                  self.other_count]
        labels = ["Exon", "CDS", "mRNA", "Other"]
        explode = (0.02, 0.01, 0.01, 0.16)
        plt.pie(values, explode=explode, labels=labels, startangle=90,
                shadow=True)
        plt.title("The number of different kinds of sequences in the "
                  "C. Elegans")
        plt.show()
        tkinter.messagebox.showinfo("Pie Diagram", "Hello!")

def main():
    headers = read_fasta_header()
    d_seq = read_fasta_seq()
    ConcencusCheck(headers, d_seq)
    exon_count = counting_exons()
    cds_count = counting_cds()
    mrna_count = counting_mrna()
    total_count = counting_total()
    other_count = calculating_other(exon_count, cds_count, mrna_count,
                                    total_count)
    GUI(exon_count, cds_count, mrna_count, total_count,
        other_count)


main()
