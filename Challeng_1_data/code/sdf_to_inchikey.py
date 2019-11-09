import rdkit
from rdkit import Chem
from rdkit.Chem import inchi
import sys
import gzip

molecules = Chem.SDMolSupplier(sys.argv[1])

csv = open(sys.argv[1] + ".inchikey", "w")

for mol in molecules:
	if mol:
		csv.write(inchi.MolToInchiKey(mol) + " " + mol.GetProp("_Name") + "\n")

csv.close()


