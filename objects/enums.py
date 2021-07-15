from enum import Enum, IntEnum, auto


class Hydropathy(Enum):
    HYDROPHOBIC = auto()
    HYDROPHILIC = auto()
    NEUTRAL = auto()


class Charge(Enum):
    POSITIVE = auto()
    NEGATIVE = auto()
    NEUTRAL = auto()


class Volume(IntEnum):
    VERY_SMALL = auto()
    SMALL = auto()
    MEDIUM = auto()
    LARGE = auto()
    VERY_LARGE = auto()


class Chemical(Enum):
    ALIPHATIC = auto()
    AROMATIC = auto()
    SULFUR = auto()
    HYDROXYL = auto()
    BASIC = auto()
    ACIDIC = auto()
    AMIDE = auto()


class HydroDonor(IntEnum):
    DONOR = auto()
    ACCEPTOR = auto()
    DONOR_ACCEPTOR = auto()
    NONE = auto()
