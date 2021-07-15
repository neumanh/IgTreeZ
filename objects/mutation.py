#!/usr/bin/python3

from objects import AminoAcid, Charge, Hydropathy, HydroDonor, Chemical, Volume


class Mutation:
    """
    Holds data on a mutation
    """
    def __init__(self, codon1: str, codon2: str, pos: int, regions: dict, no_cdr3: bool = False):
        self.codon1, self.codon2 = codon1, codon2
        self.aa1, self.aa2 = self.define_mut_by_codon(codon1, codon2)
        self.regions = regions
        self.pos = pos
        self.no_cdr3 = no_cdr3

    def __dir__(self):
        return ['pos', 'codon1', 'codon2', 'aa1', 'aa2', 'regions', 'no_cdr3']

    def __str__(self):
        # string = 'Mutation type:'
        # for att in self.__dir__():
        #     if getattr(self, att):
        #         string += f' {att}'
        string = f'Mutation: pos: {self.pos}   codon1: {self.codon1}   codon2: {self.codon2}   aa1: {self.aa1}   ' \
                 f'aa2: {self.aa2}  regions: {self.regions}  no_cdr3: {self.no_cdr3}'
        return string

    def __repr__(self):
        return f'<{str(self)}>'

    def get_one_mut_list(self):
        """
        Gets the dictionary that describe the mutation
        :return: The dictionary
        """
        arr = []  # Source
        nuc1, nuc2 = self.get_mutation_nucs(self.codon1, self.codon2)
        arr += self.define_mut_by_nuc(nuc1)

        # Trans
        arr += self.define_trans(nuc1, nuc2)

        # Region
        arr += self.define_mut_reg()

        if self.aa1 and self.aa2:
            # R/S
            if self.aa1.aa != self.aa2.aa:
                arr.append('replacement')
            else:
                arr.append('silent')

            # AA
            arr += self.define_charge()
            arr += self.define_hydro()
            arr += self.define_chemical()
            arr += self.define_polarity()
            arr += self.define_volume()
            arr += self.define_hydro_da()

        return arr

    @staticmethod
    def define_mut_by_nuc(nuc1: str):
        """
        Tests what is the source nucleotide
        :param nuc1: The first nucleotide
        :return: The mutation properties as a Mutation object
        """
        # Mutation source
        arr = []
        if nuc1 == 'A':
            arr.append('source_a')
        elif nuc1 == 'C':
            arr.append('source_c')
        elif nuc1 == 'G':
            arr.append('source_g')
        elif nuc1 == 'T':
            arr.append('source_t')
        return arr

    @staticmethod
    def define_trans(nuc1: str, nuc2: str):
        """
        Tests if the mutation is a transversion
        :param nuc1: The first nucleotide
        :param nuc2: The second nucleotide
        :return: True if it is a transversion
        """
        arr = []
        purines = ['A', 'G']
        pyrimidines = ['C', 'T']
        if ((nuc1 in purines) and (nuc2 in pyrimidines)) or ((nuc1 in pyrimidines) and (nuc2 in purines)):
            arr.append('transversion')
        else:
            arr.append('transition')

        return arr

    def define_mut_by_codon(self, codon1: str, codon2: str):
        """
        Translates the two codons
        :param codon1: The source codon
        :param codon2: The target codon
        :return: The mutation properties as a Mutation object
        """

        aa1 = self.translate_codon(codon1)
        aa2 = self.translate_codon(codon2)

        if aa1:
            aa1 = AminoAcid(aa1)
        if aa2:
            aa2 = AminoAcid(aa2)

        return aa1, aa2

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

    def define_mut_reg(self):
        """
        Defines the mutation CDR/FWR region
        :return: The updated mutation object
        """
        arr = []

        if self.pos < self.regions['cdr1']:
            arr.append('fwr1')
            arr.append('FWR')
        elif (self.pos >= self.regions['cdr1']) and (self.pos < self.regions['fwr2']):
            arr.append('cdr1')
            arr.append('CDR')
        elif (self.pos >= self.regions['fwr2']) and (self.pos < self.regions['cdr2']):
            arr.append('fwr2')
            arr.append('FWR')
        elif (self.pos >= self.regions['cdr2']) and (self.pos < self.regions['fwr3']):
            arr.append('cdr2')
            arr.append('CDR')
        elif (self.pos >= self.regions['fwr3']) and (self.pos < self.regions['cdr3']):
            arr.append('fwr3')
            arr.append('FWR')
        elif (self.pos >= self.regions['cdr3']) and (self.pos <= self.regions['cdr3_end']) and (not self.no_cdr3):
            arr.append('cdr3')
            arr.append('CDR')
        return arr

    @staticmethod
    def get_mutation_nucs(codon1: str, codon2: str):
        """
        Gets the mutated nucleotides
        :param codon1: The source codon
        :param codon2: The target codon
        :return: The two nucleotides
        """
        i = 0
        nuc1, nuc2 = None, None
        if codon1 and codon2:
            for temp_nuc in codon1:
                if temp_nuc != codon2[i]:
                    nuc1, nuc2 = temp_nuc, codon2[i]
                i += 1
        return nuc1, nuc2

    def define_charge(self):
        """
        Defines the 7 charge definitions: positive/negative/neutral for source and traget, charge change/keep
        :return: A dictionary with the definitions
        """
        arr = []
        # for aa1:
        if self.aa1.charge == Charge.POSITIVE:
            arr.append('positive_source')
        elif self.aa1.charge == Charge.NEGATIVE:
            arr.append('negative_source')
        elif self.aa1.charge == Charge.NEUTRAL:
            arr.append('neutral_charge_source')

        # for aa2:
        if self.aa2.charge == Charge.POSITIVE:
            arr.append('positive_target')
        elif self.aa2.charge == Charge.NEGATIVE:
            arr.append('negative_target')
        elif self.aa2.charge == Charge.NEUTRAL:
            arr.append('neutral_charge_target')

        # For both
        if (self.aa1.charge and self.aa2.charge) and (self.aa1.aa != self.aa2.aa):  # A replacement mutation
            if self.aa1.charge == self.aa2.charge:
                arr.append('charge_keep')
            else:
                arr.append('charge_change')

        return arr

    def define_hydro(self):
        """
        Defines the 7 hydropathy definitions: positive/negative/neutral for source and traget, charge change/keep
        :return: A dictionary with the definitions
        """
        arr = []        # for aa1:
        if self.aa1.hydro == Hydropathy.HYDROPHOBIC:
            arr.append('hydrophobic_source')
        elif self.aa1.hydro == Hydropathy.HYDROPHILIC:
            arr.append('hydrophilic_source')
        elif self.aa1.hydro == Hydropathy.NEUTRAL:
            arr.append('neutral_hydro_source')

        # for aa2:
        if self.aa2.hydro == Hydropathy.HYDROPHOBIC:
            arr.append('hydrophobic_target')
        elif self.aa2.hydro == Hydropathy.HYDROPHILIC:
            arr.append('hydrophilic_target')
        elif self.aa2.hydro == Hydropathy.NEUTRAL:
            arr.append('neutral_hydro_target')

        # For both
        if (self.aa1.hydro and self.aa2.hydro) and (self.aa1.aa != self.aa2.aa):  # A replacement mutation
            if self.aa1.hydro == self.aa2.hydro:
                arr.append('hydro_keep')
            else:
                arr.append('hydro_change')

        return arr

    def define_volume(self):
        """
        Defines the volume definitions.
        :return: A dictionary with the definitions
        """
        arr = []        # for aa1:
        if self.aa1.volume == Volume.VERY_SMALL:
            arr.append('vs_volume_source')
        elif self.aa1.volume == Volume.SMALL:
            arr.append('small_volume_source')
        elif self.aa1.volume == Volume.MEDIUM:
            arr.append('medium_volume_source')
        elif self.aa1.volume == Volume.LARGE:
            arr.append('large_volume_source')
        elif self.aa1.volume == Volume.VERY_LARGE:
            arr.append('vl_volume_source')

        # for aa2:
        if self.aa2.volume == Volume.VERY_SMALL:
            arr.append('vs_volume_target')
        elif self.aa2.volume == Volume.SMALL:
            arr.append('small_volume_target')
        elif self.aa2.volume == Volume.MEDIUM:
            arr.append('medium_volume_target')
        elif self.aa2.volume == Volume.LARGE:
            arr.append('large_volume_target')
        elif self.aa2.volume == Volume.VERY_LARGE:
            arr.append('vl_volume_target')

        # For both
        if (self.aa1.volume and self.aa2.volume) and (self.aa1.aa != self.aa2.aa):  # A replacement mutation
            if self.aa1.volume > self.aa2.volume:
                arr.append('volume_decrease')
            elif self.aa1.volume < self.aa2.volume:
                arr.append('volume_increase')
            else:
                arr.append('volume_keep')
        #
        # if not (self.aa1.volume and self.aa2.volume):
        #     print(self.codon1, '\n', self.codon2)  # TEMP

        return arr

    def define_chemical(self):
        """
        Defines the chemical definitions.
        :return: A dictionary with the definitions
        """
        arr = []        # for aa1:
        if self.aa1.chemical == Chemical.AMIDE:
            arr.append('amide_source')
        elif self.aa1.chemical == Chemical.ACIDIC:
            arr.append('acidic_source')
        elif self.aa1.chemical == Chemical.BASIC:
            arr.append('basic_source')
        elif self.aa1.chemical == Chemical.HYDROXYL:
            arr.append('hydroxyl_source')
        elif self.aa1.chemical == Chemical.SULFUR:
            arr.append('sulfur_source')
        elif self.aa1.chemical == Chemical.AROMATIC:
            arr.append('aromatic_source')
        elif self.aa1.chemical == Chemical.ALIPHATIC:
            arr.append('aliphatic_source')

        # for aa2:
        if self.aa2.chemical == Chemical.AMIDE:
            arr.append('amide_target')
        elif self.aa2.chemical == Chemical.ACIDIC:
            arr.append('acidic_target')
        elif self.aa2.chemical == Chemical.BASIC:
            arr.append('basic_target')
        elif self.aa2.chemical == Chemical.HYDROXYL:
            arr.append('hydroxyl_target')
        elif self.aa2.chemical == Chemical.SULFUR:
            arr.append('sulfur_target')
        elif self.aa2.chemical == Chemical.AROMATIC:
            arr.append('aromatic_target')
        elif self.aa2.chemical == Chemical.ALIPHATIC:
            arr.append('aliphatic_target')

        # For both
        if (self.aa1.chemical and self.aa2.chemical) and (self.aa1.aa != self.aa2.aa):  # A replacement mutation
            if self.aa1.chemical == self.aa2.chemical:
                arr.append('chemical_keep')
            else:
                arr.append('chemical_change')

        return arr

    def define_hydro_da(self):
        """
        Defines the hydrogen donor/acceptor definitions
        :return: A dictionary with the definitions
        """
        arr = []
        if self.aa1.hydro_da == HydroDonor.DONOR:  # for aa1:
            arr.append('hydro_donor_source')
        elif self.aa1.hydro_da == HydroDonor.ACCEPTOR:
            arr.append('hydro_acceptor_source')
        elif self.aa1.hydro_da == HydroDonor.DONOR_ACCEPTOR:
            arr.append('hydro_donor_acceptor_source')
        elif self.aa1.hydro_da == HydroDonor.NONE:
            arr.append('hydro_da_none_source')

        # for aa2:
        if self.aa2.hydro_da == HydroDonor.DONOR:
            arr.append('hydro_donor_target')
        elif self.aa2.hydro_da == HydroDonor.ACCEPTOR:
            arr.append('hydro_acceptor_target')
        elif self.aa2.hydro_da == HydroDonor.DONOR_ACCEPTOR:
            arr.append('hydro_donor_acceptor_target')
        elif self.aa2.hydro_da == HydroDonor.NONE:
            arr.append('hydro_da_none_target')

        # For both
        if (self.aa1.hydro_da and self.aa2.hydro_da) and (self.aa1.aa != self.aa2.aa):  # A replacement mutation
            if self.aa1.hydro_da == self.aa2.hydro_da:
                arr.append('hydro_da_keep')
            else:
                arr.append('hydro_da_change')

        return arr

    def define_polarity(self):
        """
        Defines the hydrogen donor/acceptor definitions
        :return: A dictionary with the definitions
        """
        arr = []
        if self.aa1.polarity:  # for aa1:
            arr.append('polar_source')
        else:
            arr.append('non_polar_source')

        # for aa2:
        if self.aa2.polarity:
            arr.append('polar_target')
        else:
            arr.append('non_polar_target')

        # For both
        if (self.aa1.polarity and self.aa2.polarity) and (self.aa1.aa != self.aa2.aa):  # A replacement mutation
            if self.aa1.polarity == self.aa2.polarity:
                arr.append('polarity_keep')
            else:
                arr.append('polarity_change')

        return arr
