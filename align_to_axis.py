"""Align protein to arbitrary axis with atom selection at origin.

Author: Seth Axen
E-mail: seth.axen@gmail.com
"""
from __future__ import division
import ast

from pymol import cmd, stored
from pymol.cgo import *
import numpy as np


AXES = {"x": [1., 0., 0.], "y": [0., 1., 0.], "z": [0., 0., 1.]}


def get_com(selection, state=1, mass=None, quiet=1):
    """
 DESCRIPTION

    Calculates the center of mass

    Author: Sean Law
    Michigan State University
    slaw (at) msu . edu
    """
    quiet = int(quiet)

    totmass = 0.0
    if mass is not None and not quiet:
        print "Calculating mass-weighted COM"

    state = int(state)
    model = cmd.get_model(selection, state)
    x, y, z = 0, 0, 0
    for a in model.atom:
        if (mass is not None):
            m = a.get_mass()
            x += a.coord[0] * m
            y += a.coord[1] * m
            z += a.coord[2] * m
            totmass += m
        else:
            x += a.coord[0]
            y += a.coord[1]
            z += a.coord[2]

    if (mass is not None):
        return x / totmass, y / totmass, z / totmass
    else:
        return x / len(model.atom), y / len(model.atom), z / len(model.atom)

cmd.extend("get_com", get_com)


def as_unit(v):
    """Convert to unit numpy vector"""
    v = np.asarray(v, dtype=np.float)
    return v / np.linalg.norm(v)


def create_rot_matrix(v0, v1):
    """Create 3x3 matrix of rotation from v0 onto v1 such that R*v0=v1."""
    u0 = as_unit(v0)
    u1 = as_unit(v1)
    o = np.cross(u0, u1)
    sin_ang = np.linalg.norm(o)
    cos_ang = np.dot(u0, u1)
    skew = np.array([[ 0.  , -o[2],   o[1] ],
                     [ o[2],  0.  ,  -o[0] ],
                     [-o[1],  o[0],   0.   ]], dtype=np.float)
    rot = (cos_ang * np.identity(3, dtype=np.float) + sin_ang * skew +
           (1 - cos_ang) * np.outer(o, o))
    return rot


def create_trans_matrix(rot_matrix, pre_trans_vect=np.zeros(3),
                        post_trans_vect=np.zeros(3)):
    """Create PyMOL formatted transformation matrix.

    See http://www.pymolwiki.org/index.php/Transform_selection for description
    of transformation matrix."""
    trans_matrix = np.ones((4, 4), dtype=np.float)
    trans_matrix[:3, :3] = rot_matrix
    trans_matrix[:3, 3] = post_trans_vect
    trans_matrix[3, :3] = pre_trans_vect
    return trans_matrix


def align_to_axis(obj, atom_sel=None, axis="y", axes_dict=AXES):
    """
DESCRIPTION

    Align protein to arbitrary axis with atom selection at origin.

USAGE

    align_to_axis obj [, atom_sel [, axis ]]

NOTES

    "atom_sel" must contain a single atom in "obj". If not specified, defaults
    to N-terminal nitrogen. "axis" must be either x, y, or z, or coordinates
    for an axis (e.g. [1, 2, 3]).

EXAMPLES

    # align to y-axis with N-terminus at origin
    align_to_axis obj
    # align to x-axis
    align_to_axis obj, axis=x
    # align to y-axis with C-terminus at origin
    cmd.select("cterm", "index %d:%d" % cmd.get_model("obj and not het", 1).get_residues()[-1])
    align_to_axis obj, cterm and n. c
    """
    if atom_sel is None:
        stored.list = []
        cmd.iterate("%s and not hetatm and n. n" % obj,
                    "stored.list.append(resi)")
        atom_sel = "%s and resi %s and n. n" % (obj, stored.list[0])

    if axis.lower() in axes_dict:
        target_vect = np.asarray(axes_dict[axis], dtype=np.float)
    else:
        try:
            target_vect = np.asarray(ast.literal_eval(axis), dtype=np.float)
        except (ValueError, TypeError):
            print(
                "AlignToAxisError: Provided axis is of unknown format: %s." % (
                    repr(axis)))
            return False

    cmd.reset()

    origin_coord_list = cmd.get_model(atom_sel, 1).get_coord_list()
    if len(origin_coord_list) != 1:
        print(
            "AlignToAxisError: atom selection should contain exactly 1 atom. Selection contains %d atoms." % (
                len(origin_coord_list)))
        return False
    origin_coord = np.asarray(origin_coord_list[0], dtype=np.float)

    com_coord = np.asarray(get_com(obj), dtype=float)
    com_vect = com_coord - origin_coord
    pre_trans_vect = -com_coord
    post_trans_vect = target_vect * np.linalg.norm(com_vect)
    rot_matrix = create_rot_matrix(com_vect, target_vect)
    trans_matrix = create_trans_matrix(rot_matrix, pre_trans_vect,
                                       post_trans_vect)
    cmd.transform_selection(obj, trans_matrix.flatten().tolist(),
                            homogenous=0)
    com_vect = np.asarray(get_com(obj), dtype=float)

    cmd.reset()

cmd.extend("align_to_axis", align_to_axis)
cmd.auto_arg[0]['align_to_axis'] = [cmd.object_sc, 'object', '']
cmd.auto_arg[1]['align_to_axis'] = [cmd.selection_sc, 'selection', '']
cmd.auto_arg[2]['align_to_axis'] = [lambda: cmd.Shortcut(['x', 'y', 'z']),
                                    'axis', '']
