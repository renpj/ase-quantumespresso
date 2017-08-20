# coding: utf-8
# Copyright (c) PJ Ren and Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import re
import six

from monty.io import zopen

from monty.re import regrep
from collections import defaultdict

from pymatgen.core.periodic_table import Element
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.util.io_utils import clean_lines

"""
This module implements input and output processing for PWSCF.
"""
__author__ = "PJ Ren"
__copyright__ = "Copyright 2017, PJ Ren"
__version__ = "0.1"
__maintainer__ = "PJ Ren"
__email__ = "openrpj@gmail.com"
__date__ = "8/19/2017"

def get_atom_mass(symbol):
    n = ase.data.atomic_numbers[symbol]
    return ase.data.atomic_masses[n-1]

class PWInput(object):
    """
    Base input file class. Right now, only supports no symmetry and is
    very basic. Initially adopted from pymatgen
    """

    def __init__(self, atoms, pseudo=None, control=None, system=None,
                 electrons=None, ions=None, cell=None, kpoints_mode="automatic",
                 kpoints_grid=(1, 1, 1),kpoints_shift=(0, 0, 0)):
        """
        Initializes a PWSCF input file.

        Args:
            atoms (ase.Atoms): Input atoms.
            pseudo (dict): A dict of the pseudopotentials to use. Default to None.
            control (dict): Control parameters. Refer to official PWSCF doc
                on supported parameters. Default to {"calculation": "scf"}
            system (dict): System parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            electrons (dict): Electron parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            ions (dict): Ions parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            cell (dict): Cell parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            kpoints_mode (str): Kpoints generation mode. Default to automatic.
            kpoints_grid (sequence): The kpoint grid. Default to (1, 1, 1).
            kpoints_shift (sequence): The shift for the kpoints. Defaults to
                (0, 0, 0).
        """
        self.atoms = atoms
        sections = {}
        sections["control"] = control or {"calculation": "scf"}
        sections["system"] = system or {}
        sections["electrons"] = electrons or {}
        sections["ions"] = ions or {}
        sections["cell"] = cell or {}
        if pseudo == None:
            raise PWInputError("Missing all pseudo specification!")
        else:
            for symbol in set(self.atoms.get_chemical_symbols()):
                if not in pseudo:
                    raise PWInputError("Missing %s in pseudo specification!" 
                                       % symbol)
        self.pseudo = pseudo
        self.sections = sections
        self.kpoints_mode = kpoints_mode
        self.kpoints_grid = kpoints_grid
        self.kpoints_shift = kpoints_shift

    def __str__(self):
        out = []

        def to_str(v):
            if isinstance(v, six.string_types):
                return "'%s'" % v
            elif isinstance(v, float):
                return "%s" % str(v).replace("e", "d")
            elif isinstance(v, bool):
                if v:
                    return ".TRUE."
                else:
                    return ".FALSE."
            return v

        for k1 in ["control", "system", "electrons", "ions", "cell"]:
            v1 = self.sections[k1]
            out.append("&%s" % k1.upper())
            sub = []
            for k2 in sorted(v1.keys()):
                if isinstance(v1[k2], list): # for Hubbard_U etc.
                    n = 1
                    for l in v1[k2][:len(self.pseudo)]: 
                        sub.append("  %s(%d) = %s" % (k2, n, to_str(v1[k2][n-1])))
                        n += 1
                else:
                    sub.append("  %s = %s" % (k2, to_str(v1[k2])))
            if k1 == "system":
                if 'ibrav' not in self.sections[k1]:
                    sub.append("  ibrav = 0")
                if 'nat' not in self.sections[k1]:
                    sub.append("  nat = %d" % len(self.atoms))
                if 'ntyp' not in self.sections[k1]:
                    sub.append("  ntyp = %d" % len(self.pseudo))
            sub.append("/")
            out.append(",\n".join(sub))

        out.append("ATOMIC_SPECIES")
        for k, v in sorted(self.pseudo.items(), key=lambda i: i[0]):
            e = re.match(r"[A-Z][a-z]?", k).group(0)
            if self.pseudo is not None:
                p = v
            else:
                p = v['pseudo']
            out.append("  %s  %.4f %s" % (k, get_atom_mass(e), p)) 

        out.append("ATOMIC_POSITIONS crystal")
        for name,pos in zip(self.atoms.get_chemical_symbols(),self.atoms.get_scaled_positions()):
            out.append("  %s %.6f %.6f %.6f" % (name, pos[0],pos[1],pos[2]))

        out.append("K_POINTS %s" % self.kpoints_mode)
        kpt_str = ["%s" % i for i in self.kpoints_grid]
        kpt_str.extend(["%s" % i for i in self.kpoints_shift])
        out.append("  %s" % " ".join(kpt_str))
        out.append("CELL_PARAMETERS angstrom")
        for vec in self.atoms.cell:
            out.append("  %f %f %f" % (vec[0], vec[1], vec[2]))
        return "\n".join(out)

    def write_file(self, filename):
        """
        Write the PWSCF input file.

        Args:
            filename (str): The string filename to output to.
        """
        with open(filename, "w") as f:
            f.write(self.__str__())

    @staticmethod
    def from_file(filename):
        """
        Reads an PWInput object from a file.

        Args:
            filename (str): Filename for file

        Returns:
            PWInput object
        """
        with zopen(filename, "rt") as f:
            return PWInput.from_string(f.read())

    @staticmethod
    def from_string(string):
        """
        Reads an PWInput object from a string.

        Args:
            string (str): PWInput string

        Returns:
            PWInput object
        """
        lines = list(clean_lines(string.splitlines()))

        def input_mode(line):
            if line[0] == "&":
                return ("sections", line[1:].lower())
            elif "ATOMIC_SPECIES" in line:
                return ("pseudo", )
            elif "K_POINTS" in line:
                return ("kpoints", line.split("{")[1][:-1])
            elif "CELL_PARAMETERS" in line or "ATOMIC_POSITIONS" in line:
                return ("structure", line.split("{")[1][:-1])
            elif line == "/":
                return None
            else:
                return mode

        sections = {"control": {}, "system": {}, "electrons": {}, 
                    "ions": {}, "cell":{}}
        pseudo = {}
        pseudo_index = 0
        lattice = []
        species = []
        coords = []
        structure = None
        site_properties = {"pseudo":[]}
        mode = None
        for line in lines:
            mode = input_mode(line)
            if mode == None:
                pass
            elif mode[0] == "sections":
                section = mode[1]
                m = re.match(r'(\w+)\(?(\d*?)\)?\s*=\s*(.*)', line)
                if m:
                    key = m.group(1).strip()
                    key_ = m.group(2).strip()
                    val = m.group(3).strip()
                    if key_ != "":
                        if sections[section].get(key, None) == None:
                            val_ = [0.0]*20 # MAX NTYP DEFINITION
                            val_[int(key_)-1] = PWInput.proc_val(key, val)
                            sections[section][key] = val_

                            site_properties[key] = []
                        else:
                            sections[section][key][int(key_)-1] = PWInput.proc_val(key, val) 
                    else:
                        sections[section][key] = PWInput.proc_val(key, val)

            elif mode[0] == "pseudo":
                m = re.match(r'(\w+)\s+(\d*.\d*)\s+(.*)', line)
                if m:
                    pseudo[m.group(1).strip()] = {}
                    pseudo[m.group(1).strip()]["index"] = pseudo_index
                    pseudo[m.group(1).strip()]["pseudopot"] = m.group(3).strip()
                    pseudo_index += 1
            elif mode[0] == "kpoints":
                m = re.match(r'(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
                if m:
                    kpoints_grid = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
                    kpoints_shift = (int(m.group(4)), int(m.group(5)), int(m.group(6)))
                else:
                    kpoints_mode = mode[1]
            elif mode[0] == "structure":
                m_l = re.match(r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)', line)
                m_p = re.match(r'(\w+)\s+(-?\d+\.\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)', line)
                if m_l:
                    lattice += [ float(m_l.group(1)), float(m_l.group(2)), float(m_l.group(3)) ]
                elif m_p:
                    site_properties["pseudo"].append(pseudo[m_p.group(1)]["pseudopot"])
                    species += [pseudo[m_p.group(1)]["pseudopot"].split(".")[0]]
                    coords += [[float(m_p.group(2)), float(m_p.group(3)), float(m_p.group(4))]]

                    for k, v in site_properties.items():
                        if k != "pseudo":
                            site_properties[k].append(sections['system'][k][pseudo[m_p.group(1)]["index"]])
                if mode[1] == "angstrom":
                    coords_are_cartesian = True
                elif mode[1] == "crystal":
                    coords_are_cartesian = False

        if coords_are_cartesian:
            atoms = ase.Atoms(symbols=species,positions=coords,info=site_properties)
        else:
            atoms = ase.Atoms(symbols=species,scaled_positions=coords,info=site_properties)
        return PWInput(atoms=atoms, control=sections["control"],
                       system=sections["system"], electrons=sections["electrons"], 
                       ions=sections["ions"], cell=sections["cell"], kpoints_mode=kpoints_mode,
                       kpoints_grid=kpoints_grid, kpoints_shift=kpoints_shift)

    def proc_val(key, val):
        """
        Static helper method to convert PWINPUT parameters to proper type, e.g.,
        integers, floats, etc.

        Args:
            key: PWINPUT parameter key
            val: Actual value of PWINPUT parameter.
        """
        float_keys = ('etot_conv_thr','forc_conv_thr','conv_thr','Hubbard_U','Hubbard_J0','defauss',
                      'starting_magnetization',)

        int_keys = ('nstep','iprint','nberrycyc','gdir','nppstr','ibrav','nat','ntyp','nbnd','nr1',
                    'nr2','nr3','nr1s','nr2s','nr3s','nspin','nqx1','nqx2','nqx3','lda_plus_u_kind',
                    'edir','report','esm_nfit','space_group','origin_choice','electron_maxstep',
                    'mixing_ndim','mixing_fixed_ns','ortho_para','diago_cg_maxiter','diago_david_ndim',
                    'nraise','bfgs_ndim','if_pos','nks','nk1','nk2','nk3','sk1','sk2','sk3','nconstr')

        bool_keys = ('wf_collect','tstress','tprnfor','lkpoint_dir','tefield','dipfield','lelfield',
                     'lorbm','lberry','lfcpopt','monopole','nosym','nosym_evc','noinv','no_t_rev',
                     'force_symmorphic','use_all_frac','one_atom_occupations','starting_spin_angle',
                     'noncolin','x_gamma_extrapolation','lda_plus_u','lspinorb','london',
                     'ts_vdw_isolated','xdm','uniqueb','rhombohedral','realxz','block',
                     'scf_must_converge','adaptive_thr','diago_full_acc','tqr','remove_rigid_rot',
                     'refold_pos')

        def smart_int_or_float(numstr):
            if numstr.find(".") != -1 or numstr.lower().find("e") != -1:
                return float(numstr)
            else:
                return int(numstr)

        try:
            if key in bool_keys:
                if val.lower() == ".true.":
                    return True
                elif val.lower() == ".false.":
                    return False
                else:
                    raise ValueError(key + " should be a boolean type!")

            if key in float_keys:
                return float(re.search(r"^-?\d*\.?\d*d?-?\d*", val.lower()).group(0).replace("d", "e"))

            if key in int_keys:
                return int(re.match(r"^-?[0-9]+", val).group(0))

        except ValueError:
            pass

        try:
            val = val.replace("d","e")
            return smart_int_or_float(val)
        except ValueError:
            pass

        if "true" in val.lower():
            return True
        if "false" in val.lower():
            return False

        m = re.match(r"^[\"|'](.+)[\"|']$", val)
        if m:
            return m.group(1)



class PWInputError(BaseException):
    pass


class PWOutput(object):

    patterns = {
        "energies": r'total energy\s+=\s+([\d\.\-]+)\sRy',
        "ecut": r'kinetic\-energy cutoff\s+=\s+([\d\.\-]+)\s+Ry',
        "lattice_type": r'bravais\-lattice index\s+=\s+(\d+)',
        "celldm1": r"celldm\(1\)=\s+([\d\.]+)\s",
        "celldm2": r"celldm\(2\)=\s+([\d\.]+)\s",
        "celldm3": r"celldm\(3\)=\s+([\d\.]+)\s",
        "celldm4": r"celldm\(4\)=\s+([\d\.]+)\s",
        "celldm5": r"celldm\(5\)=\s+([\d\.]+)\s",
        "celldm6": r"celldm\(6\)=\s+([\d\.]+)\s",
        "nkpts": r"number of k points=\s+([\d]+)"
    }

    def __init__(self, filename):
        self.filename = filename
        self.data = defaultdict(list)
        self.read_pattern(PWOutput.patterns)
        for k, v in self.data.items():
            if k == "energies":
                self.data[k] = [float(i[0][0]) for i in v]
            elif k in ["lattice_type", "nkpts"]:
                self.data[k] = int(v[0][0][0])
            else:
                self.data[k] = float(v[0][0][0])

    def read_pattern(self, patterns, reverse=False,
                     terminate_on_match=False, postprocess=str):
        """
        General pattern reading. Uses monty's regrep method. Takes the same
        arguments.

        Args:
            patterns (dict): A dict of patterns, e.g.,
                {"energy": r"energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)"}.
            reverse (bool): Read files in reverse. Defaults to false. Useful for
                large files, esp OUTCARs, especially when used with
                terminate_on_match.
            terminate_on_match (bool): Whether to terminate when there is at
                least one match in each key in pattern.
            postprocess (callable): A post processing function to convert all
                matches. Defaults to str, i.e., no change.

        Renders accessible:
            Any attribute in patterns. For example,
            {"energy": r"energy\\(sigma->0\\)\\s+=\s+([\\d\\-.]+)"} will set the
            value of self.data["energy"] = [[-1234], [-3453], ...], to the
            results from regex and postprocess. Note that the returned
            values are lists of lists, because you can grep multiple
            items on one line.
        """
        matches = regrep(self.filename, patterns, reverse=reverse,
                         terminate_on_match=terminate_on_match,
                         postprocess=postprocess)
        self.data.update(matches)

    def get_celldm(self, i):
        return self.data["celldm%d" % i]

    @property
    def final_energy(self):
        return self.data["energies"][-1]

    @property
    def lattice_type(self):
        return self.data["lattice_type"]

