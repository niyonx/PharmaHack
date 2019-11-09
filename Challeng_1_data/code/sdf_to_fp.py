import sys
import rdkit
from rdkit import Chem
from rdkit.Chem.rdmolops import *

molecules = Chem.SDMolSupplier(sys.argv[1])

for mol in molecules:
	if mol:
		#fp = RDKFingerprint(mol, fpSize=2048)
		fp = RDKFingerprint(mol, fpSize=2048, minPath=1, maxPath=7, nBitsPerHash=2, useHs=1)
		print(str(fp.ToBase64()) + "\t" + str(mol.GetProp("_Name")))


