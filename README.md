# pymol_scripts
Scripts for performing useful tasks within [PyMOL](https://www.pymol.org/ "PyMOL")

##Setup##

Assuming this repo is located in your home directory, add the following to your `~/.pymolrc` file:

```python
python

import os
import glob

script_repo_path = os.path.join(os.path.expanduser("~"), "pymol_scripts")
scripts = glob.glob(script_repo_path+"/*.py")
for script in scripts:
    cmd.run(script)

python end
```

When PyMOL starts up, all of these scripts should be available in PyMOL.

##Script Categories##

###Protein Alignment###
- `align_to_axis`
    + Transform protein so that an atom selection (default N-terminal nitrogen) is located at the origin and the center of mass is located along a provided axis (default y-axis).

###Small Molecule Alignment###
_Note: These scripts require that the RDKit distribution be in your `$PYTHONPATH` and compiled for the same version of Python as PyMOL._
- `rdkitalign`
    + Align molecules using [RDKit](rdkit.org "RDKit")'s 3D molecular alignment. [See documentation](http://www.rdkit.org/Python_Docs/rdkit.Chem.rdMolAlign-module.html#AlignMol "RDKit AlignMol Documentation").
- `open3dalign`
    + Align molecules using [RDKit](rdkit.org "RDKit")'s interface to [Open3DAlign](http://open3dalign.sourceforge.net/ "Open3DAlign"). [See documentation](http://www.rdkit.org/Python_Docs/rdkit.Chem.rdMolAlign-module.html#GetO3A "RDKit GetO3A Documentation").
