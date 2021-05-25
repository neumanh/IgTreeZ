#!/usr/bin/python3


class AminoAcid:
    """
    Holds data on a mutation
    """
    def __init__(self):
        self.positive = False
        self.negative = False
        self.hydrophilic = False
        self.hydrophobic = False

    def __dir__(self):
        return ['positive', 'negative', 'hydrophilic', 'hydrophobic']

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

