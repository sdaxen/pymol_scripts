"""Align two molecules using RDKit.
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


def rdkitalign(mobile, target, mobile_state=-1, target_state=-1,
               reflect=0, max_iters=500):
    """
DESCRIPTION

    Align two molecules using RDKit.

USAGE

    rdkitalign mobile, target, [, target_state [, mobile_state [, reflect [,
        max_iters ]]]]

NOTES

    "rdkitalign" uses RDKit's alignment interface. "target_state"
    and "mobile_state" indicate the states of the molecules that will be used
    the alignment, though all states will be transformed accordingly. These
    default to -1, the use of the currently visible states.
    "reflect" reflects the conformation of the mobile molecule (defaults to 0).

EXAMPLES

    # align molA onto molB
    rdkitalign molA, molB
    """
    print("Aligning %s to %s with RDKit" % (mobile, target))
    mobile_mol = obj_to_rdmol(mobile, int(mobile_state))
    target_mol = obj_to_rdmol(target, int(target_state))
    try:
        rmsd, trans_matrix = rdMolAlign.GetAlignmentTransform(
            mobile_mol, target_mol, reflect=bool(int(reflect)),
            maxIters=int(max_iters))
    except RuntimeError:
        print " RDKitAlign-Error: alignment failed."
        return False
    cmd.transform_selection(mobile, trans_matrix.flatten().tolist(),
                            homogenous=1)
    if HAS_NP:
        inv_array_string = "["+", ".join(
            ["%.4g" % x for x
             in np.linalg.inv(trans_matrix).flatten().tolist()])+"]"
        print("Molecules are aligned. To reverse the alignment, run:")
        print('cmd.transform_selection("%s", %s, homogenous=1)' % (
            mobile, inv_array_string))

    print("RMSD: %.4f" % rmsd)

cmd.extend("rdkitalign", rdkitalign)
cmd.auto_arg[0]['rdkitalign'] = [cmd.object_sc, 'object', '']
cmd.auto_arg[1]['rdkitalign'] = [cmd.object_sc, 'object', '']
