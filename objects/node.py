#!/usr/bin/python3


class Node:
    """
    Holds data on transitions
    """
    def __init__(self, name: str, dist: int):
        self.name = name
        self.dist = dist # Reperesent the distance to the father

    def __str__(self):
        string = f'name = {self.name} distance = {self.dist}'
        return string