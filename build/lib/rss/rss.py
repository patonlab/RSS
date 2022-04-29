#!/usr/bin/env python

###########################################################################################.
###########################################################################################
###                                                                                     ###
###            Script to calculate RSS for Quantum Chemistry output files               ###
###                                                                                     ###
###  Authors: Shree Sowndarya S. V.,                                                    ###
###  Please, report any bugs or suggestions to:                                         ###
###  svss@colostate.edu                                                                 ###
###                                                                                     ###
###########################################################################################
###########################################################################################.

import numpy as np
import pandas as pd
import os, sys, argparse
import dbstep.Dbstep as db
import glob
from openbabel import openbabel as ob
import cclib

from rdkit import Chem

periodictable = [
    "",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Uub",
    "Uut",
    "Uuq",
    "Uup",
    "Uuh",
    "Uus",
    "Uuo",
]


def output_to_spin_bv(files, name_csv, type):
    all_data = pd.DataFrame()
    counter = 0

    for log in files:

        name = log.split(".")[0]

        data = cclib.io.ccread(log)
        spins = data.atomspins[type]
        atoms = data.atomnos

        i = 1
        for spin, atom in zip(spins, atoms):
            all_data.at[counter, "file"] = log
            all_data.at[counter, "name"] = name
            all_data.at[counter, "atom_sym"] = periodictable[atom]
            all_data.at[counter, "atom_idx"] = i

            all_data.at[counter, "spin_density"] = spin

            if periodictable[atom] != "H":
                sterics_scan = db.dbstep(log, atom1=i, volume=True)
                all_data.at[counter, "buried_vol"] = float(sterics_scan.bur_vol)
            else:
                all_data.at[counter, "buried_vol"] = None

            counter += 1
            i += 1

    all_data = all_data[all_data.atom_sym != "H"]
    all_data["fractional_spin"] = all_data.groupby("name")["spin_density"].apply(
        lambda x: x.abs() / x.abs().sum()
    )

    max_cdf = all_data.loc[all_data.groupby("name")["fractional_spin"].idxmax()]
    max_cdf["stability"] = max_cdf["buried_vol"] + 50 * (1 - max_cdf["fractional_spin"])

    all_data.to_csv(name_csv + "_all.csv", index=False)
    max_cdf.to_csv(name_csv + ".csv", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--files",
        help="Computational output files to calculate buried volume and spin",
        nargs="+",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output name for CSVs of spin and buried volume",
        default="output",
    )
    parser.add_argument(
        "-t",
        "--type",
        help="Type of spin to use (Default: Mulliken spin)",
        default="mulliken",
        choices=["mulliken", "lowdin"],
    )
    args = parser.parse_args()

    output_to_spin_bv(args.files, args.output, args.type)


if __name__ == "__main__":
    main()
