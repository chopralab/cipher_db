from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

smiles_list = [
    "O=C(N[C@@H](C)CC1=C(/C=C/C(OCC)=O)SC=C1)NC[C@@H](N(C)C)CC2=CC=C(O)C=C2",
    "O=C(N[C@@H](C)CC1=C(/C=C/C(CC)=O)SC=C1)NC[C@@H](N(C)C)CC2=CC=C(O)C=C2",
    "O=C(N[C@@H](C)CC1=CSC=C1)NC[C@@H](N(C)C)CC2=CC(/C=C/S(F)(=O)=O)=C(O)C=C2",
    "O=C(N[C@@H](C)CC1=C(/C=C/C2=CC=CC=C2)SC=C1)NC[C@@H](N(C)C)CC3=CC=C(O)C=C3",
    "O=C(N[C@@H](C)CC1=C(/C=C/C2=CC=C(Cl)C=C2)SC=C1)NC[C@@H](N(C)C)CC3=CC=C(O)C=C3",
    "O=C(N[C@@H](C)CC1=CSC=C1)NC[C@@H]([NH+](C)C)CC2=CC=C(OS(=O)(F)=O)C=C2",
    "O=C(N[C@@H](C)CC1=CSC=C1)NC[C@@H](N(C)C)CC2=CC=C(OS(NCCCC)(=O)=O)C=C2",
    "O=C(N[C@@H](C)CC1=CSC=C1)NC[C@@H](N(C)C)CC2=CC=C(OS(N(C)CC3=CC=CC=C3)(=O)=O)C=C2",
    "O=C(N[C@@H](C)CC1=CSC=C1)NC[C@@H](N(C)C)CC2=CC=C(OS(NC3CCCCC3)(=O)=O)C=C2"
]

mol_form_list = []

for smi in smiles_list:
    mol = Chem.MolFromSmiles(smi)
    mol_form_list.append(CalcMolFormula(mol))

print(mol_form_list)