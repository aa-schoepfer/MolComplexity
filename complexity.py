import argparse
import numpy as np
import networkx as nx
from itertools import combinations
from typing import Dict, TextIO, List, Tuple
from numpy.typing import ArrayLike

atomic_symbols: Dict[int, str] = {
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "Ne",
    11: "Na",
    12: "Mg",
    13: "Al",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    18: "Ar",
    19: "K",
    20: "Ca",
    21: "Sc",
    22: "Ti",
    23: "V",
    24: "Cr",
    25: "Mn",
    26: "Fe",
    27: "Co",
    28: "Ni",
    29: "Cu",
    30: "Zn",
    31: "Ga",
    32: "Ge",
    33: "As",
    34: "Se",
    35: "Br",
    36: "Kr",
    37: "Rb",
    38: "Sr",
    39: "Y",
    40: "Zr",
    41: "Nb",
    42: "Mo",
    43: "Tc",
    44: "Ru",
    45: "Rh",
    46: "Pd",
    47: "Ag",
    48: "Cd",
    49: "In",
    50: "Sn",
    51: "Sb",
    52: "Te",
    53: "I",
    54: "Xe",
    55: "Cs",
    56: "Ba",
    57: "La",
    58: "Ce",
    59: "Pr",
    60: "Nd",
    61: "Pm",
    62: "Sm",
    63: "Eu",
    64: "Gd",
    65: "Tb",
    66: "Dy",
    67: "Ho",
    68: "Er",
    69: "Tm",
    70: "Yb",
    71: "Lu",
    72: "Hf",
    73: "Ta",
    74: "W",
    75: "Re",
    76: "Os",
    77: "Ir",
    78: "Pt",
    79: "Au",
    80: "Hg",
    81: "Tl",
    82: "Pb",
    83: "Bi",
    84: "Po",
    85: "At",
    86: "Rn",
    87: "Fr",
    88: "Ra",
    89: "Ac",
    90: "Th",
    91: "Pa",
    92: "U",
    93: "Np",
    94: "Pu",
    95: "Am",
    96: "Cm",
    97: "Bk",
    98: "Cf",
    99: "Es",
    100: "Fm",
    101: "Md",
    102: "No",
    103: "Lr",
    104: "Rf",
    105: "Db",
    106: "Sg",
    107: "Bh",
    108: "Hs",
    109: "Mt",
    110: "Ds",
    111: "Rg",
    112: "Cn",
    113: "Nh",
    114: "Fl",
    115: "Mc",
    116: "Lv",
    117: "Ts",
    118: "Og",
}

atomic_numbers: Dict[str, int] = {v: k for k, v in atomic_symbols.items()}


class MolReader:
    def __init__(self, stream: TextIO, strict: bool = True) -> None:
        self.strict = strict
        self.g = self._get_graph(stream)

    def _get_stereo(self, G, bond_orders, bond_stereos):
        pass

    def _get_graph(self, stream: TextIO) -> None:
        lines: List[str] = stream.read().splitlines()

        version: str = " V2000"
        if len(lines[3]) != 39 or lines[3][33:39] != version:
            raise ValueError(
                f"Molfile V2000 expected. Counts line length of 39 expcted, got {len(lines[3])}. '{version}' expcted, got {lines[3][33:39]}"
            )

        end_tag: str = "M  END"
        if lines[-1] != end_tag:
            raise ValueError(f"No '{end_tag}' end tag found. Got {lines[-1]} instead.")

        n_atoms: int = int(lines[3][0:3])
        n_bonds: int = int(lines[3][3:6])

        atoms: ArrayLike[int] = np.array(
            [atomic_numbers[line[31:34].strip()] for line in lines[4 : 4 + n_atoms]]
        )
        bonds: List[Tuple[int, int]] = [
            (int(line[0:3]) - 1, int(line[3:6]) - 1)
            for line in lines[4 + n_atoms : 4 + n_atoms + n_bonds]
        ]

        bond_orders: ArrayLike[int] = np.array(
            [int(line[6:9]) for line in lines[4 + n_atoms : 4 + n_atoms + n_bonds]]
        )

        bond_stereos: ArrayLike[int] = np.array(
            [int(line[9:12]) for line in lines[4 + n_atoms : 4 + n_atoms + n_bonds]]
        )

        G = nx.Graph()
        G.add_edges_from(bonds)

        # TODO
        # atom_stereos: ArrayLike[int] = self._get_stereo(G, bond_orders, bond_stereos)

        for i, n in enumerate(G.nodes):
            G.nodes[n]["an"] = atoms[i]

        for i, (s, e) in enumerate(G.edges):
            G[s][e]["bo"] = bond_orders[i]

        for i, (s, e) in enumerate(G.edges):
            G[s][e]["ch"] = bond_stereos[i]

        return G


def complexity(G) -> float:
    assert G.number_of_nodes() > 1, "Single atoms don't work."
    nmax = -1
    bmax = -1
    for i in reversed(range(1, G.number_of_edges())):
        Gs = []
        for j in combinations(G.edges, i):
            g = G.edge_subgraph(j)
            if nx.is_connected(g):
                Gs.append(g)
        uniques = []
        for g in Gs:
            wl1 = nx.weisfeiler_lehman_graph_hash(g, node_attr="an", edge_attr="bo")
            wl2 = nx.weisfeiler_lehman_graph_hash(g, edge_attr="ch")
            if wl1 + wl2 not in uniques:
                uniques.append(wl1 + wl2)
        lu = len(uniques)
        # print(i, lu)
        if lu < nmax:
            break
        if lu >= nmax:
            nmax = lu
            bmax = i
    # print("###")
    # print(f"bmax = {bmax}, nmax = {nmax}")
    if bmax == 1:
        return 0
    else:
        return np.log(nmax) / np.log(bmax)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="MolComplexity",
        description="Calculate molecular complexity from Mol a file.",
        epilog='"Fools ignore complexity. Pragmatists suffer it. Some can avoid it. Geniuses remove it." --Alan Perlis',
    )

    parser.add_argument("filename", help="A V2000 Molfile")
    # parser.add_argument('-f', '--flexible', default=False, help="Allow some flexibility in Molfiles sythax. ")
    flexible = False
    args = parser.parse_args()

    with open(args.filename) as file:
        mol = MolReader(file, strict=flexible)

    print(complexity(mol.g))
