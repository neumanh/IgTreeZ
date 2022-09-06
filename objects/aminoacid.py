#!/usr/bin/python3

from objects import Charge, Hydropathy, HydroDonor, Chemical, Volume


class AminoAcid:
    """
    Holds data on a mutation
    """

    def __init__(self, aa: str):
        self.aa = aa
        self.charge = self.define_charge()
        self.hydro = self.define_hydropathy()
        self.volume = self.define_volume()
        self.chemical = self.define_chemical()
        self.hydro_da = self.define_hydro_donor()
        self.polarity = self.define_polar()

    def __dir__(self):
        return ['aa', 'charge', 'hydro', 'volume', 'chemical', 'hydro_da', 'polarity']

    def __str__(self):
        # string = f'Amino acid type: {self.aa}'
        # for att in self.__dir__():
        #     if getattr(self, att):
        #         string += f' {att}'
        string = f'AA: {self.aa}   charge: {self.charge}   hydro: {self.hydro}   volume: {self.volume}   ' \
                 f'chemical: {self.chemical} hydro_da: {self.hydro_da} polarity: {self.polarity}'
        return string

    def __repr__(self):
        return f'<{str(self)}>'

    def get_dict(self):
        """
        Creates a dictionary of the mutation
        :return: The dictionary, in which each attrubutes is saved as key
        """
        dic = {}
        for att in self.__dir__():
            dic[att] = getattr(self, att)

        return dic

    def define_hydropathy(self):
        """
        Define the Hydropathy of the amino-acid
        :return:
        """
        # By IMGT's definition: http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
        hydrophilic_aa = ['R', 'N', 'D', 'E', 'Q', 'K']
        hydrophobic_aa = ['A', 'V', 'I', 'L', 'F', 'W', 'M', 'C']

        hydro = None
        if self.aa and (self.aa != '_'):  # Not a None or a stop codon
            if self.aa in hydrophilic_aa:
                hydro = Hydropathy.HYDROPHILIC
            elif self.aa in hydrophobic_aa:
                hydro = Hydropathy.HYDROPHOBIC
            else:
                hydro = Hydropathy.NEUTRAL
        return hydro

    def define_charge(self):
        """
        Define the charge of the amino-acid
        :return:
        """
        # By IMGT's definition: http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
        positive_aa = ['R', 'K', 'H']
        negative_aa = ['D', 'E']

        charge = None
        if self.aa and (self.aa != '_'):  # Not a None or a stop codon
            if self.aa in positive_aa:
                charge = Charge.POSITIVE
            elif self.aa in negative_aa:
                charge = Charge.NEGATIVE
            else:
                charge = Charge.NEUTRAL

        return charge

    def define_volume(self):
        """
        Define the Volume of the amino-acid
        :return:
        """
        # By IMGT's definition: http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
        vol_vs = ['A', 'G', 'S']
        vol_s = ['N', 'D', 'C', 'P', 'T']
        vol_m = ['Q', 'E', 'H', 'V']
        vol_l = ['R', 'I', 'L', 'K', 'M']
        vol_vl = ['F', 'W', 'Y']

        vol = None

        if self.aa and (self.aa != '_'):  # Not a None or a stop codon
            if self.aa in vol_vs:
                vol = Volume.VERY_SMALL
            elif self.aa in vol_s:
                vol = Volume.SMALL
            elif self.aa in vol_m:
                vol = Volume.MEDIUM
            elif self.aa in vol_l:
                vol = Volume.LARGE
            elif self.aa in vol_vl:
                vol = Volume.VERY_LARGE
        return vol

    def define_chemical(self):
        """
        Define the Chemical of the amino-acid
        :return:
        """
        # By IMGT's definition: http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
        chem_aliphatic = ['A', 'G', 'I', 'L', 'P', 'V']
        chem_aromatic = ['F', 'W', 'Y']
        chem_sulfur = ['C', 'M']
        chem_hydroxyl = ['S', 'T']
        chem_basic = ['R', 'H', 'K']
        chem_acidic = ['D', 'E']
        chem_amide = ['N', 'Q']

        chem = None
        if self.aa and (self.aa != '_'):  # Not a None or a stop codon
            if self.aa in chem_aliphatic:
                chem = Chemical.ALIPHATIC
            elif self.aa in chem_aromatic:
                chem = Chemical.AROMATIC
            elif self.aa in chem_sulfur:
                chem = Chemical.SULFUR
            elif self.aa in chem_hydroxyl:
                chem = Chemical.HYDROXYL
            elif self.aa in chem_basic:
                chem = Chemical.BASIC
            elif self.aa in chem_acidic:
                chem = Chemical.ACIDIC
            elif self.aa in chem_amide:
                chem = Chemical.AMIDE

        return chem

    def define_hydro_donor(self):
        """
        Define the charge of the amino-acid
        :return:
        """
        # By IMGT's definition: http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
        hydro_donor = ['R', 'K', 'W']
        hydro_accept = ['D', 'E']
        hydro_donacc = ['N', 'Q', 'H', 'S', 'T', 'Y']
        hydro_none = ['A', 'C', 'G', 'I', 'L', 'M', 'F', 'P', 'V']

        hd = None
        if self.aa and (self.aa != '_'):  # Not a None or a stop codon
            if self.aa in hydro_donor:
                hd = HydroDonor.DONOR
            elif self.aa in hydro_accept:
                hd = HydroDonor.ACCEPTOR
            elif self.aa in hydro_donacc:
                hd = HydroDonor.DONOR_ACCEPTOR
            elif self.aa in hydro_none:
                hd = HydroDonor.NONE
        return hd

    def define_polar(self):
        """
        Define the polarity of the amino-acid
        :return:
        """
        # By IMGT's definition: http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
        # nonpolar_aa = ['A', 'C', 'G', 'I', 'L', 'M', 'F', 'P', 'W', 'V']

        polar_aa = ['R', 'N', 'D', 'Q', 'E', 'H', 'K', 'S', 'T', 'Y']

        polarity = None
        if self.aa and (self.aa != '_'):  # Not a None or a stop codon
            if self.aa in polar_aa:
                polarity = True
            else:
                polarity = False
        return polarity
