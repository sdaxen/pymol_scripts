"""Align two molecules using Open3DAlign.
Author: Seth Axen
E-mail: seth.axen@gmail.com
"""
import os
import tempfile
from pymol import cmd
import rdkit
import rdkit.Chem
from rdkit.Chem import rdMolAlign
try:
    import numpy as np
    HAS_NP = True
except ImportError:
    HAS_NP = False


def obj_to_rdmol(obj, state=-1):
    mol_file = tempfile.mktemp(suffix='.mol', prefix=obj)
    cmd.save(mol_file, obj, state)
    mol = rdkit.Chem.MolFromMolFile(mol_file)
    os.remove(mol_file)
    return mol


def open3dalign(mobile, target, mobile_state=-1, target_state=-1,
                options=0, max_iters=500):
    """
DESCRIPTION

    Align two molecules using Open3DAlign.

USAGE

    open3dalign mobile, target, [, target_state [, mobile_state [, options [,
        max_iters ]]]]

NOTES

    "open3dalign" uses RDKit's interface for Open3DAlign. "target_state"
    and "mobile_state" indicate the states of the molecules that will be used
    the alignment, though all states will be transformed accordingly. These
    default to -1, the use of the currently visible states.
    For description of other parameters, see RDKit's Open3DAlign documentation:
    http://www.rdkit.org/Python_Docs/rdkit.Chem.rdMolAlign-module.html#GetO3A

EXAMPLES

    # align molA onto molB
    open3dalign molA, molB
    """
    print("Aligning %s to %s with Open3DAlign" % (mobile, target))
    mobile_mol = obj_to_rdmol(mobile, int(mobile_state))
    target_mol = obj_to_rdmol(target, int(target_state))
    alignment = rdMolAlign.GetO3A(mobile_mol, target_mol, options=int(options),
                                  maxIters=int(max_iters))
    rmsd, trans_matrix = alignment.Trans()
    cmd.transform_selection(mobile, trans_matrix.flatten().tolist(),
                            homogenous=1)
    if HAS_NP:
        inv_array_string = "["+", ".join(
            ["%.4g" % x for x
             in np.linalg.inv(trans_matrix).flatten().tolist()])+"]"
        print("Molecules are aligned. To reverse the alignment, run:")
        print('cmd.transform_selection("%s", %s, homogenous=1)' % (
            mobile, inv_array_string))

    score = alignment.Score()
    print("RMSD: %.4f" % rmsd)
    print("Score: %.4f" % score)

cmd.extend("open3dalign", open3dalign)
cmd.auto_arg[0]['open3dalign'] = [cmd.object_sc, 'object', '']
cmd.auto_arg[1]['open3dalign'] = [cmd.object_sc, 'object', '']
