# coding: utf-8
# Copyright (c) PJ Ren.
# Distributed under the terms of the MIT License.
import ase

def get_atom_mass(symbol):
    n = ase.data.atomic_numbers[symbol]
    return ase.data.atomic_masses[n-1]

def clean_lines(string_list, remove_empty_lines=True):
    """
    Strips whitespace, carriage returns and empty lines from a list of strings.

    Args:
        string_list: List of strings
        remove_empty_lines: Set to True to skip lines which are empty after
            stripping.

    Returns:
        List of clean strings with no whitespaces.
    """

    for s in string_list:
        clean_s = s
        if '#' in s:
            ind = s.index('#')
            clean_s = s[:ind]
        clean_s = clean_s.strip()
        if (not remove_empty_lines) or clean_s != '':
            yield clean_s
