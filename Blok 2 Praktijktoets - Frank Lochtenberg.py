# Made by Frank Lochtenberg
# Version 8

import re
import matplotlib.pyplot as plt
import tkinter
from tkinter import messagebox


def read_fasta():
    """ Reads the fasta file and puts the lines from the file in a list.

        :return: list of lines from the fasta file.
        """
    try:
        fasta = open("Caenorhabditis_elegans.cds_pep.all.fa", "r")
        fasta_list = []
        for line in fasta:
            fasta_list.append(line)
        return fasta_list
    except FileNotFoundError:
        print("The fasta file was not found in read_fasta().")


def read_gff3():
    """ Reads the gff3 file and puts the lines from the file in a list.

    :return: list of lines from the gff3 file.
    """
    try:
        gff3 = open("Caenorhabditis_elegans.gff3", "r")
        gff3_list = []
        for line in gff3:
            gff3_list.append(line)
        return gff3_list
    except FileNotFoundError:
        print("The gff3 file was not found in read_gff3().")


def fasta_header(fasta_list):
    """ Puts the amino acids headers from the fasta file in a list.

    :param fasta_list: list of lines from the fasta file.
    :return: list of headers of amino acids.
    """
    try:
        headers = []
        for line in fasta_list:
            if ">" in line:
                headers.append(line)
        return headers
    except SyntaxError:
        print("There was a Syntax Error in fasta_header()")


def fasta_seq(fasta_list):
    """ Puts the amino acid sequences from the fasta file in a list.

    :param fasta_list: list of lines from the fasta file.
    :return: list of amino acids sequences.
    """
    try:
        t_seq = []
        d_seq = []
        for line in fasta_list:
            line = line.strip()
            if ">" in line:
                if t_seq != []:
                    d_seq.append("".join(t_seq))
                    t_seq = []
            else:
                t_seq.append(line)
        d_seq.append("".join(t_seq))
        return d_seq
    except SyntaxError:
        print("There was a Syntax Error in fasta_seq()")


def counting_exons(gff3_list):
    """Counts the frequency of exons from the gff3 list.

    :param gff3_list: list of lines from the gff3 file.
    :return: An integer which is the frequency of exons in the gff3 file
    """
    try:
        exon_count = 0
        for line in gff3_list:
            if "x" in line[11:14]:
                exon_count += 1
        return exon_count
    except IndexError:
        print("There was an Index Error in counting_exons()")


def counting_cds(gff3_list):
    """Counts the frequency CDS from the gff3 list.

    :param gff3_list: list of lines from the gff3 file.
    :return: An integer which is the frequency of CDS in the gff3 file.
    """
    try:
        cds_count = 0
        for line in gff3_list:
            if "D" in line[11:14]:
                cds_count += 1
        return cds_count
    except IndexError:
        print("There was an Index Error in counting_cds()")


def counting_mrna(gff3_list):
    """Counts the frequency of mRNA from the gff3 list.

    :param gff3_list: list of lines from the gff3 file.
    :return: An integer which is the frequency of mRNA in the gff3 file.
    """
    try:
        mrna_count = 0
        for line in gff3_list:
            if "m" in line[11:14]:
                mrna_count += 1
        return mrna_count
    except IndexError:
        print("There was an Index Error in counting_mrna()")


def counting_total(gff3_list):
    """Counts the total lines from the gff3 list.

    :param gff3_list: list of lines from the gff3 file.
    :return: An integer which is the total of lines in the gff3 file.
    """
    try:
        total_count = 0
        for line in gff3_list:
            if line != "":
                total_count += 1
        return total_count
    except IndexError:
        print("There was an Attribute Error in counting_total()")


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
    try:
        other_count = total_count - exon_count - cds_count - mrna_count
        return other_count
    except TypeError:
        print("There was a Type Error in calculating_other()")


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
        try:
            self.set_headers(headers)
            self.set_aa_seq(seq)
            self.set_consensus()
        except AttributeError:
            print("There was an Attribute Error in Class ConsensusCheck"
                  " in the init()")

    def set_headers(self, headers):
        """ Converting the headers list into an object.

        :param headers: list of headers of amino acids.
        :return: nothing
        """
        try:
            self.headers = headers
        except AttributeError:
            print("There was an Attribute Error in Class ConsensusCheck"
                  " set_headers()")

    def get_headers(self):
        """ Returns the headers list object.

        :return: An headers list object
        """
        try:
            return self.headers
        except AttributeError:
            print("There was an Attribute Error in Class ConsensusCheck"
                  " in get_headers()")

    def set_aa_seq(self, seq):
        """ Converting the amino acid sequences list into an
        object.

        :param seq: an amino acid sequences list object
        :return: nothing
        """
        try:
            self.aa_seq = seq
        except AttributeError:
            print("There was an Attribute Error in Class ConsensusCheck"
                  " in set_aa_seq()")

    def get_aa_seq(self):
        """ Returns the amino acid sequences list object.

        :return:
        """
        try:
            return self.aa_seq
        except AttributeError:
            print("There was an Attribute Error in Class ConsensusCheck"
                  " in get_aa_seq()")

    def set_consensus(self):
        """ Searches for the zinc finger consensus and puts the
        sequences which contain the consensus in a list and then into a
        object.

        :return: nothing
        """
        try:
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
            self.ncss = ncss
        except AttributeError:
            print("There was an Attribute Error in Class ConsensusCheck"
                  " in set_consensus()")
        except ModuleNotFoundError:
            print("The re-module could not be found in Class "
                  "ConsensusCheck in set_consensus()")

    def get_consensus(self):
        """ Returns the zinc finger consensus list object

        :return: an zinc finger consensus list object
        """
        try:
            return self.consensus
        except AttributeError:
            print("There was an Attribute Error in Class ConsensusCheck"
                  " in get_consensus()")

    def get_non_consensus(self):
        """ Returns an integer object which is the number of amino acid
        sequences which do not contain the zinc finger consensus.

        :return: An integer object which is the number of amino acid
        sequences which do not contain the zinc finger consensus.
        """
        try:
            return self.ncss
        except AttributeError:
            print("There was an Attribute Error in Class ConsensusCheck"
                  " in get_non_consensus()")


class GUI:
    def __init__(self, exon_count, cds_count, mrna_count, total_count,
                 other_count):
        try:
            # Converting the given parameters into a object
            self.exons_count = exon_count
            self.cds_count = cds_count
            self.mrna_count = mrna_count
            self.total_count = total_count
            self.other_count = other_count

            # Transferring the zinc finger consensus list object and an
            # integer object which is the number of amino acid sequences
            # which do not contain the zinc finger consensus to this
            # class and making it a new variable and putting it in a new
            # object
            fasta_list = read_fasta()
            headers = fasta_header(fasta_list)
            seq = fasta_seq(fasta_list)
            zfcs = ConsensusCheck(headers, seq).get_consensus()
            ncss = ConsensusCheck(headers, seq).get_non_consensus()
            self.ncss = ncss
            self.zfcs = zfcs

            # Making a integer object of the number of zinc finger
            # consensus there are
            zinc_finger_count = 0
            for _ in self.zfcs:
                zinc_finger_count += 1
            self.zfct = zinc_finger_count

            # Makes the main window
            self.main_window = tkinter.Tk()
            self.main_window.title("Information over the sequences")

            # Makes frames
            self.intro_frame = tkinter.Frame(self.main_window)
            self.info_buttons_frame = tkinter.Frame(self.main_window)
            self.graph_buttons_frame = tkinter.Frame(self.main_window)
            self.css_buttons_frame = tkinter.Frame(self.main_window)
            self.hello_button_frame = tkinter.Frame(self.main_window)
            self.quit_frame = tkinter.Frame(self.main_window)

            # Makes a label for the introduction sentence
            self.intro_label = tkinter.Label(self.intro_frame,
                                             text="Click the buttons "
                                                  "for information!")

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
                                               text="Total of items "
                                                    "Info",
                                               command=self.total)

            # Places the buttons for info
            self.exon_button.pack(side="left")
            self.cds_button.pack(side="left")
            self.mrna_button.pack(side="left")
            self.other_button.pack(side="left")
            self.total_button.pack(side="left")

            # Makes the graph buttons
            self.seq_pie_button = tkinter.Button(self.
                                                 graph_buttons_frame,
                                                 text="Show Sequences "
                                                      "Pie Diagram",
                                                 command=self.
                                                 pie_diagram_seq)
            self.bar_button = tkinter.Button(self.graph_buttons_frame,
                                             text="Show Bar Diagram",
                                             command=self.bar_diagram)
            self.stacked_bar_buton = tkinter.Button(self.
                                                    graph_buttons_frame,
                                                    text="Show "
                                                         "Stacked "
                                                         "Bar "
                                                         "Diagram",
                                                    command=self.
                                                    stacked_bar_graph)

            # Places the graph buttons
            self.seq_pie_button.pack(side="left")
            self.bar_button.pack(side="left")
            self.stacked_bar_buton.pack(side="left")

            # Makes the consensus buttons
            self.consensus_button = tkinter.Button(self.
                                                   css_buttons_frame,
                                                   text="Zinc Finger "
                                                        "Consensus "
                                                        "Frequency "
                                                        "Info",
                                                   command=self.
                                                   consensus)
            self.css_pie = tkinter.Button(self.css_buttons_frame,
                                          text="Show Consensus Pie "
                                               "Diagram",
                                          command=self.pie_diagram_css)

            # Places the consensus buttons
            self.consensus_button.pack(side="left")
            self.css_pie.pack(side="left")

            # Make Hello Button
            self.hello_button = tkinter.Button(self.hello_button_frame,
                                               text="Hello",
                                               command=self.hello)

            # Place Hello Button
            self.hello_button.pack(side="left")

            # Makes the quit button
            self.quit_button = tkinter.Button(self.quit_frame,
                                              text="Quit",
                                              command=self.main_window.
                                              destroy)

            # Places the quit button
            self.quit_button.pack(side="right")

            # Places the frames in the GUI.
            self.intro_frame.pack()
            self.hello_button_frame.pack()
            self.info_buttons_frame.pack()
            self.graph_buttons_frame.pack()
            self.css_buttons_frame.pack()
            self.quit_frame.pack()

            # Makes the main window visible.
            tkinter.mainloop()
        except TypeError:
            print("There was a Type Error in Class GUI in the init()")
        except SyntaxError:
            print("There was a Syntax Error in Class GUI "
                  "in the init()")
        except NameError:
            print("There was a Name Error in Class GUI in the init()")
        except AttributeError:
            print("There was an Attribute Error in Class GUI "
                  "in the init()")
        except ModuleNotFoundError:
            print("The tkinter-module could not be found in Class GUI "
                  "in the init()")

    def exons(self):
        """ Shows the information over the exons sequences in a
        messagebox.

        :return: nothing
        """
        try:
            tkinter.messagebox.showinfo("Exons Info",
                                        "There are " +
                                        str(self.exons_count) +
                                        " exon sequences in the "
                                        "C. Elegans.")
        except AttributeError:
            print("There was an Attribute Error in Class GUI "
                  "in exons()")
        except ModuleNotFoundError:
            print("The tkinter messagebox-module could not be found in "
                  "Class GUI in the exons()")

    def cds(self):
        """ Shows the information over the CDS sequences in a
        messagebox.

        :return: nothing
        """
        try:
            tkinter.messagebox.showinfo("CDS Info",
                                        "There are " +
                                        str(self.cds_count) +
                                        " CDS sequences in the C. "
                                        "Elegans.")
        except AttributeError:
            print("There was an Attribute Error in Class GUI "
                  "in cds()")
        except ModuleNotFoundError:
            print("The tkinter messagebox-module could not be found in "
                  "Class GUI in the cds()")

    def mrna(self):
        """ Shows the information over the mRNA sequences in a
        messagebox.

        :return: nothing
        """
        try:
            tkinter.messagebox.showinfo("mRNA Info",
                                        "There are " +
                                        str(self.mrna_count) +
                                        " mRNA sequences in the "
                                        "C. Elegans.")
        except AttributeError:
            print("There was an Attribute Error in Class GUI "
                  "in mrna()")
        except ModuleNotFoundError:
            print("The tkinter messagebox-module could not be found in "
                  "Class GUI in the mrna()")

    def other(self):
        """ Shows the information over the other items than exon, CDS &
        mRNA sequences in a messagebox.

        :return: nothing
        """
        try:
            tkinter.messagebox.showinfo("Other items Info",
                                        "There are " +
                                        str(self.other_count) +
                                        " of other items than exon, CDS"
                                        " & mRNA sequences in the "
                                        "C. Elegans.")
        except AttributeError:
            print("There was an Attribute Error in Class GUI "
                  "in other()")
        except ModuleNotFoundError:
            print("The tkinter messagebox-module could not be found in "
                  "Class GUI in the other()")

    def total(self):
        """ Shows the information over all sequences in a messagebox.

        :return: nothing
        """
        try:
            tkinter.messagebox.showinfo("Total of items Info",
                                        "There are " +
                                        str(self.total_count) +
                                        " of total sequences in the "
                                        "C. Elegans.")
        except AttributeError:
            print("There was an Attribute Error in Class GUI "
                  "in the total()")
        except ModuleNotFoundError:
            print("The tkinter messagebox-module could not be found in "
                  "Class GUI in the total()")

    def consensus(self):
        """ Shows the information over the zinc finger consensus in a
        messagebox.

        :return: nothing
        """
        try:
            tkinter.messagebox.showinfo("Zinc Finger Consensus Info",
                                        "There are " +
                                        str(self.zfct) +
                                        " of the Zinc Finger Consensus "
                                        "constated in the C. Elegans.")
        except AttributeError:
            print("There was an Attribute Error in Class GUI "
                  "in consensus()")
        except ModuleNotFoundError:
            print("The tkinter messagebox-module could not be found in "
                  "Class GUI in the consensus()")

    def pie_diagram_seq(self):
        """ Makes a pie diagram over all the sequences.

        :return: nothing
        """
        try:
            values = [self.exons_count, self.cds_count, self.mrna_count,
                      self.other_count]
            labels = ["Exon", "CDS", "mRNA", "Other"]
            explode = (0.02, 0.01, 0.01, 0.16)
            plt.pie(values, explode=explode, labels=labels,
                    startangle=90, shadow=True,
                    colors=["y", "r", "b", "k"])
            plt.title("The number of different kinds of sequences in "
                      "the C. Elegans")
            plt.show()
        except AttributeError:
            print("There was an Attribute Error in Class GUI "
                  "in pie_diagram_set()")
        except ModuleNotFoundError:
            print("The matplotlib-module could not be found in "
                  "Class GUI in the pie_diagram_seq()")

    def bar_diagram(self):
        """ Makes a bar diagram over all the sequences.

        :return: nothing
        """
        try:
            x = ["Exon", "CDS", "mRNA", "Other"]
            y = [self.exons_count, self.cds_count, self.mrna_count,
                 self.other_count]
            plt.bar(x, y, color=["y", "r", "b", "k"])
            plt.title("The number of different kinds of sequences in "
                      "the C. Elegans")
            plt.xlabel("Kind of sequences")
            plt.ylabel("Frequency of kind of sequences")
            plt.grid(True, color="w", axis="y")
            plt.show()
        except AttributeError:
            print("There was an Attribute Error in Class GUI "
                  "in bar_diagram()")
        except ModuleNotFoundError:
            print("The matplotlib-module could not be found in "
                  "Class GUI in the bar_diagram()")

    def stacked_bar_graph(self):
        """ Makes a stacked bar diagram over all the sequences

        :return: nothing
        """
        try:
            x = ["Exons, CDS & mRNA sequences", "Other sequences"]
            e = [self.exons_count, 0]
            c = [self.cds_count, 0]
            m = [self.mrna_count, 0]
            o = [0, self.other_count]
            plt.bar(x, m, color="b")
            plt.bar(x, c, color="r", bottom=m)
            plt.bar(x, e, color="y", bottom=c)
            plt.bar(x, o, color="k")
            plt.title("The number of different kinds of sequences in "
                      "the C. Elegans")
            plt.xlabel("All Sequences")
            plt.ylabel("The frequency of different kinds of sequences "
                       "in the C. Elegans")
            legend = ["The frequency of exon sequences",
                      "The frequency of CDS sequences",
                      "The frequency of mRNA sequences",
                      "The frequency of other kinds of sequences"]
            plt.legend(legend)
            plt.grid(True, color="w", axis="y")
            plt.show()
        except AttributeError:
            print("There was an Attribute Error in Class GUI "
                  "in stacked_bar_graph()")
        except ModuleNotFoundError:
            print("The matplotlib-module could not be found in "
                  "Class GUI in the stacked_bar_graph()")

    def pie_diagram_css(self):
        """ Makes a pie diagram over the difference between sequences
        with the zinc finger consensus and without.

        :return: nothing
        """
        try:
            values = [self.zfct, self.ncss]
            labels = ["Zinc finger consensus sequences frequency",
                      "Non zinc finger consensus sequences frequency"]
            plt.pie(values, labels=labels, startangle=90, shadow=True,
                    colors=["y", "k"])
            plt.title("The difference between sequences with the zinc "
                      "finger consensus found in the sequences and "
                      "without the consensus")
            plt.show()
        except AttributeError:
            print("There was a Attribute Error in Class GUI "
                  "in pie_diagram_css()")
        except ModuleNotFoundError:
            print("The matplotlib-module could not be found in "
                  "Class GUI in the pie_diagram_css()")

    def hello(self):
        """ Shows hello in the messagebox.

        :return: Sadly no hello back...
        """
        try:
            tkinter.messagebox.showinfo(":)", "I just wanted to say "
                                              "Hello!")
        except ModuleNotFoundError:
            print("The tkinter messagebox-module could not be found in "
                  "Class GUI in the total()")


def main():
    fasta_list = read_fasta()
    gff3_list = read_gff3()
    headers = fasta_header(fasta_list)
    d_seq = fasta_seq(fasta_list)
    ConsensusCheck(headers, d_seq)
    exon_count = counting_exons(gff3_list)
    cds_count = counting_cds(gff3_list)
    mrna_count = counting_mrna(gff3_list)
    total_count = counting_total(gff3_list)
    other_count = calculating_other(exon_count, cds_count, mrna_count,
                                    total_count)
    GUI(exon_count, cds_count, mrna_count, total_count,
        other_count)


main()
