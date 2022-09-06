#!/usr/bin/python3


class Treat:
    """
    Holds lineage tree attributes
    """
    def __init__(self):
        self.id = 0
        self.nodes = 0
        self.leaves = 0
        self.od_avg = 0
        self.sod_avg = 0
        self.rootd = 0
        self.drsn_min = 0
        self.dasn_min = 0
        self.trunk = 0
        self.dlfsn_avg = 0
        self.pl_min = 0


    def __dir__(self):
        return ['id', 'nodes', 'leaves', 'od_avg', 'sod_avg', 'rootd', 'drsn_min', 'dasn_min', 'trunk', 'dlfsn_avg', 'pl_min']

    def __str__(self):
        string = 'Tree attributes:'
        for att in self.__dir__():
            if getattr(self, att):
                string += f'\t{att}={getattr(self, att)}'

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

