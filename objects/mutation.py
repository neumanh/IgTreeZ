#!/usr/bin/python3


class Mutation:
    """
    Holds data on a mutation
    """

    def __init__(self):
        self.transition = False
        self.transversion = False
        self.r = False
        self.s = False
        self.cdr = False
        self.fwr = False
        self.source_a = False
        self.source_c = False
        self.source_g = False
        self.source_t = False
        self.source_positive = False
        self.source_negative = False
        self.dest_positive = False
        self.dest_negative = False
        self.source_hydrophilic = False
        self.source_hydrophobic = False
        self.dest_hydrophilic = False
        self.dest_hydrophobic = False
        self.cdr1 = False
        self.cdr2 = False
        self.cdr3 = False
        self.fwr1 = False
        self.fwr2 = False
        self.fwr3 = False

    def __dir__(self):
        return ['transition', 'transversion', 'r', 's', 'cdr1', 'cdr2', 'cdr3', 'fwr1', 'fwr2', 'fwr3',
                 'source_a', 'source_c', 'source_g', 'source_t',
                'source_positive', 'source_negative', 'dest_positive', 'dest_negative', 'source_hydrophilic',
                'source_hydrophobic', 'dest_hydrophilic', 'dest_hydrophobic']

    def __str__(self):
        string = 'Mutation type:'
        for att in self.__dir__():
            if getattr(self, att):
                string += f' {att}'

        return string

    def __repr__(self):
        return f'<{str(self)}>'

    def get_dict(self):
        """
        Creates a dictionary of the mutation
        :return: The dictionary, in wich each attrubutes is saved as key
        """
        dic = {}
        for att in self.__dir__():
            dic[att] = getattr(self, att)

        return dic

    def describe_mut_by_nuc(self, nuc1: str, nuc2: str):
        """
        Tests if the 2 input codons code for the different amino acid and other parameters
        :param nuc1: The first nucleotide
        :param nuc2: The second nucleotide
        :return: The mutation properties as a Mutation object
        """

        # Mutaion source
        if nuc1 == 'A':
            self.source_a = True
        elif nuc1 == 'C':
            self.source_c = True
        elif nuc1 == 'G':
            self.source_g = True
        elif nuc1 == 'T':
            self.source_t = True

        # Tansition / transversion check
        if self.is_transversion(nuc1, nuc2):
            self.transversion = True
        else:
            self.transition = True

    @staticmethod
    def is_transversion(nuc1: str, nuc2: str):
        """
        Tests if the mutation is a transversion
        :param nuc1: The first nucleotide
        :param nuc2: The second nucleotide
        :return: True if it is a transversion
        """
        purines = ['A', 'G']
        pyrimidines = ['C', 'T']

        transv = False

        if (nuc1 in purines) and (nuc2 in pyrimidines):
            transv = True
        elif (nuc1 in pyrimidines) and (nuc2 in purines):
            transv = True

        return transv

    def describe_mut_by_codon(self, codon1: str, codon2: str):
        """
        Tests if the 2 input codons code for the different amino acid and other parameters
        :param codon1: The first codon
        :param codon2: The second codon
        :return: The mutation properties as a Mutation object
        """

        aa1 = self.translate_codon(codon1)
        aa2 = self.translate_codon(codon2)

        # hydrophilic_aa = ['R', 'N', 'D', 'E', 'Q', 'K', 'S', 'T']
        # hydrophobic_aa = ['A', 'V', 'G,' 'I', 'L', 'F', 'P', 'W', 'Y']
        # positive_aa = ['R', 'K', 'H']
        # negative_aa = ['D', 'E']

        # By IMGT's definition: http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
        hydrophilic_aa = ['R', 'N', 'D', 'E', 'Q', 'K']
        hydrophobic_aa = ['A', 'V', 'I', 'L', 'F', 'W', 'M', 'C']
        positive_aa = ['R', 'K', 'H']
        negative_aa = ['D', 'E']

        if aa1 and aa2:

            # R/S
            if aa1 != aa2:
                self.r = True
            else:
                self.s = True

            # hydrophilic_aa / hydrophobic_aa
            if aa1 in hydrophilic_aa:
                self.source_hydrophilic = True
            if aa1 in hydrophobic_aa:
                self.source_hydrophobic = True
            if aa2 in hydrophilic_aa:
                self.dest_hydrophilic = True
            if aa2 in hydrophobic_aa:
                self.dest_hydrophobic = True

            # positive / negative
            if aa1 in positive_aa:
                self.source_positive = True
            if aa1 in negative_aa:
                self.source_negative = True
            if aa2 in positive_aa:
                self.dest_positive = True
            if aa2 in negative_aa:
                self.dest_negative = True

    @staticmethod
    def translate_codon(codon):
        """
        Translates 3 nucleotides to amino acid.
        From https://www.geeksforgeeks.org/dna-protein-python-3/
        :param codon: nucleotide sequence
        :return: amino acid sequence
        """

        table = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
            'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
        }
        if codon in table:
            aa = table[codon]
        else:
            aa = None
        return aa

    @staticmethod
    def get_region(pos: int, regions: dict):
        """
        Gets the region of the input position (CDR or FWR)
        :type pos: int
        :param pos: The position in the Ig sequence
        :param regions: A dictionary of the FWRs and CDRs starting positions
        :return:  'CDR' if cdr, 'FWR' if framework, or None if neither
        """
        if (pos < regions['cdr1']) or ((pos >= regions['fwr2']) and (pos < regions['cdr2'])) or \
                ((pos >= regions['fwr3']) and (pos < regions['cdr3'])):
            region = 'FWR'
        elif ((pos >= regions['cdr1']) and (pos < regions['fwr2'])) or \
                ((pos >= regions['cdr2']) and (pos < regions['fwr3'])) or \
                ((pos >= regions['cdr3']) and (pos <= regions['cdr3_end'])):
            region = 'CDR'
        else:
            region = None

        return region

    def define_mut_pos(self, pos: int, regions: dict, no_cdr3: bool = False):
        """
        Define the mutation CDR/FWR region
        :param pos: The mutation position (in the aligned sequence)
        :param regions: A dictionary of the FWRs and CDRs starting positions
        :param no_cdr3: Avoid updating the cdr3 region
        :return: The updated mutation object
        """
        if pos < regions['cdr1']:
            self.fwr1 = True
        elif (pos >= regions['cdr1']) and (pos < regions['fwr2']):
            self.cdr1 = True
        elif (pos >= regions['fwr2']) and (pos < regions['cdr2']):
            self.fwr2 = True
        elif (pos >= regions['cdr2']) and (pos < regions['fwr3']):
            self.cdr2 = True
        elif (pos >= regions['fwr3']) and (pos < regions['cdr3']):
            self.fwr3 = True
        elif (pos >= regions['cdr3']) and (pos <= regions['cdr3_end']) and (not no_cdr3):
            self.cdr3 = True
