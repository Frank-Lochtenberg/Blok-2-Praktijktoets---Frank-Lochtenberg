# Made by Frank Lochtenberg
# Version 6

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


class ConsensusCheck:
    """Checks if the is a consensus pattern in the amino acids sequences

    """

    def __init__(self, headers, seq):
        """ Converting the parameters into an object.
        Checking with the given parameters if there is a concencus
        pattern in the amino acids sequences.

        :param headers: list of headers of amino acids.
        :param seq: list of amino acids sequences.
        """
        self.set_headers(headers)
        self.set_aa_seq(seq)
        self.set_consensus()

    def set_headers(self, headers):
        """ Converting the headers list into an object.

        :param headers: list of headers of amino acids.
        :return: nothing
        """
        self.headers = headers

    def get_headers(self):
        """ Returns the headers list object.

        :return: An headers list object
        """
        return self.headers

    def set_aa_seq(self, seq):
        """ Converting the amino acid sequences list into an
        object.

        :param seq: an amino acid sequences list object
        :return: nothing
        """
        self.aa_seq = seq

    def get_aa_seq(self):
        """ Returns the amino acid sequences list object.

        :return:
        """
        return self.aa_seq

    def set_consensus(self):
        """ Searches for the zinc finger consensus and puts the
        sequences which contain the consensus in a list and then into a
        object.

        :return: a integer which is the number of amino acid sequences
        which do not contain the zinc finger consensus.
        """
        ncss = 0
        cssl = []
        for aa_seq in self.get_aa_seq():
            css = re.search(r"C[A-Z]{2}C[A-Z]{2}C[A-Z]{5}C[A-Z]{2}"
                            r"C[A-Z]{2}C", aa_seq)
            if css:
                cssl.append(css.group())
                self.consensus = cssl
            else:
                ncss += 1
        return ncss

    def get_consensus(self):
        """ Returns the zinc finger consensus list object

        :return: an zinc finger consensus list object
        """
        return self.consensus


class GUI:
    def __init__(self, exon_count, cds_count, mrna_count, total_count,
                 other_count):

        # Converting the given parameters into a object
        self.exons_count = exon_count
        self.cds_count = cds_count
        self.mrna_count = mrna_count
        self.total_count = total_count
        self.other_count = other_count

        # Transferring the zinc finger consensus list object to this
        # class and making it a new variable and putting it in a new
        # object
        headers = read_fasta_header()
        seq = read_fasta_seq()
        zfc = ConsensusCheck(headers, seq).get_consensus()
        self.zfc = zfc

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
        self.consensus_button = tkinter.Button(self.info_buttons_frame,
                                               text="Zinc Finger "
                                                    "Consensus Info",
                                               command=self.consensus)

        # Places the buttons for info
        self.exon_button.pack(side="left")
        self.cds_button.pack(side="left")
        self.mrna_button.pack(side="left")
        self.other_button.pack(side="left")
        self.total_button.pack(side="left")
        self.consensus_button.pack(side="left")

        # Makes the graph buttons
        self.pie_button = tkinter.Button(self.graph_buttons_frame,
                                         text="Show Pie Diagram",
                                         command=self.pie_diagram)
        self.bar_button = tkinter.Button(self.graph_buttons_frame,
                                         text="Show Bar Diagram",
                                         command=self.bar_diagram)
        self.stacked_bar_button = tkinter.Button(self.
                                                 graph_buttons_frame,
                                                 text="Show Stacked Bar"
                                                      " Diagram",
                                                 command=self.
                                                 stacked_bar_diagram)

        # Make Hello Button
        self.hello_button = tkinter.Button(self.graph_buttons_frame,
                                           text="Hello",
                                           command=self.hello)

        # Place Hello Button
        self.hello_button.pack(side="left")

        # Places the graph buttons
        self.pie_button.pack(side="left")
        self.bar_button.pack(side="left")
        self.stacked_bar_button.pack(side="left")

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
        """ Shows the information over the exons sequences in a
        messagebox.

        :return: nothing
        """
        tkinter.messagebox.showinfo("Exons Info",
                                    "There are " +
                                    str(self.exons_count) +
                                    " exon sequences in the "
                                    "C. Elegans.")

    def cds(self):
        """ Shows the information over the CDS sequences in a
        messagebox.

        :return: nothing
        """
        tkinter.messagebox.showinfo("CDS Info",
                                    "There are " +
                                    str(self.cds_count) +
                                    " CDS sequences in the C. Elegans.")

    def mrna(self):
        """ Shows the information over the mRNA sequences in a
        messagebox.

        :return: nothing
        """
        tkinter.messagebox.showinfo("mRNA Info",
                                    "There are " +
                                    str(self.mrna_count) +
                                    " mRNA sequences in the "
                                    "C. Elegans.")

    def other(self):
        """ Shows the information over the other items than exon, CDS &
        mRNA sequences in a messagebox.

        :return: nothing
        """
        tkinter.messagebox.showinfo("Other items Info",
                                    "There are " +
                                    str(self.other_count) +
                                    " of other items than exon, CDS & "
                                    "mRNA sequences in the C. Elegans.")

    def total(self):
        """ Shows the information over all sequences in a messagebox.

        :return: nothing
        """
        tkinter.messagebox.showinfo("Total of items Info",
                                    "There are " +
                                    str(self.total_count) +
                                    " of total sequences in the "
                                    "C. Elegans.")

    def consensus(self):
        """ Shows the information over the zinc finger consensus in a
        messagebox.

        :return: nothing
        """
        zinc_finger_count = 0
        for _ in self.zfc:
            zinc_finger_count += 1
        tkinter.messagebox.showinfo("Zinc Finger Consensus Info",
                                    "There are " +
                                    str(zinc_finger_count) +
                                    " of the Zinc Finger Consensus "
                                    "constated in the C. Elegans.")

    def pie_diagram(self):
        """ Makes a pie diagram over all the sequences.

        :return: nothing
        """
        values = [self.exons_count, self.cds_count, self.mrna_count,
                  self.other_count]
        labels = ["Exon", "CDS", "mRNA", "Other"]
        explode = (0.02, 0.01, 0.01, 0.16)
        plt.pie(values, explode=explode, labels=labels, startangle=90,
                shadow=True, colors="yrbk")
        plt.title("The number of different kinds of sequences in the "
                  "C. Elegans")
        plt.show()

    def bar_diagram(self):
        """ Makes a bar diagram over all the sequences.

        :return: nothing
        """
        x = ["Exon", "CDS", "mRNA", "Other"]
        y = [self.exons_count, self.cds_count, self.mrna_count,
             self.other_count]
        plt.bar(x, y, color="yrbk")
        plt.title("The number of different kinds of sequences in the "
                  "C. Elegans")
        plt.xlabel("Kind of sequences")
        plt.ylabel("Frequency of kind of sequences")
        plt.grid(True, color="w", axis="y")
        plt.show()

    def stacked_bar_diagram(self):
        """ Makes a stacked bar diagram over all the sequences

        :return: nothing
        """
        x = ["Exons, CDS & mRNA sequences", "Other sequences"]
        e = [self.exons_count, 0]
        c = [self.cds_count, 0]
        m = [self.mrna_count, 0]
        o = [0, self.other_count]
        plt.bar(x, m, color="b")
        plt.bar(x, c, color="r", bottom=m)
        plt.bar(x, e, color="y", bottom=c)
        plt.bar(x, o, color="k")
        plt.title("The number of different kinds of sequences in the "
                  "C. Elegans")
        plt.xlabel("All Sequences", color="g")
        plt.ylabel("The frequency of different kinds of sequences in "
                   "the C. Elegans", color="g")
        legend = ["The frequency of exon sequences",
                  "The frequency of CDS sequences",
                  "The frequency of mRNA sequences",
                  "The frequency of other kinds of sequences"]
        plt.legend(legend)
        plt.grid(True, color="w", axis="y")
        plt.show()

    def hello(self):
        """ Shows hello in the messagebox.

        :return: Sadly no hello back...
        """
        tkinter.messagebox.showinfo(":)", "I just wanted to say Hello!")


def main():
    headers = read_fasta_header()
    d_seq = read_fasta_seq()
    ConsensusCheck(headers, d_seq)
    exon_count = counting_exons()
    cds_count = counting_cds()
    mrna_count = counting_mrna()
    total_count = counting_total()
    other_count = calculating_other(exon_count, cds_count, mrna_count,
                                    total_count)
    GUI(exon_count, cds_count, mrna_count, total_count,
        other_count)


main()
