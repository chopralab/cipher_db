"""adapated from https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/sascorer.py"""
import math
import sys
from typing import List

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


class SAScorer:
    def __init__(self, xss: List[List[float]]):

        self.fscores = {x: xs[0] for xs in xss for x in xs[1:]}

    def __call__(self, smi: str) -> float:
        return self.score(smi)

    def score(self, smi: str) -> float:
        mol = Chem.MolFromSmiles(smi)
        fp = rdMolDescriptors.GetMorganFingerprint(mol, 2)
        fps = fp.GetNonzeroElements()

        score1 = 0.0
        nf = 0
        for bid, v in fps.items():
            nf += v
            score1 += v * self.fscores.get(bid, -4)
        score1 /= nf

        n_atom = mol.GetNumAtoms()
        n_chiral_center = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        n_bridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
        n_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        n_macrocycle = sum(len(x) > 8 for x in mol.GetRingInfo().AtomRings())

        p_size = n_atom**1.005 - n_atom
        p_stereo = math.log10(n_chiral_center + 1)
        p_bridge = math.log10(n_bridgehead + 1)
        p_spiro = math.log10(n_spiro + 1)
        p_macrocycle = math.log10(2) if n_macrocycle > 0 else 0

        score2 = -(p_size + p_stereo + p_spiro + p_bridge + p_macrocycle)
        score3 = math.log(float(n_atom) / len(fps)) * 0.5 if n_atom > len(fps) else 0.0

        sascore = score1 + score2 + score3

        LB = -4.0
        UB = 2.5
        sascore = 11.0 - (sascore - LB + 1) / (UB - LB) * 9.0

        if sascore > 8.0:
            sascore = 8.0 + math.log(sascore + 1.0 - 9.0)
        if sascore > 10.0:
            sascore = 10.0
        elif sascore < 1.0:
            sascore = 1.0

        return sascore


if __name__ == "__main__":
    scorer = SAScorer(sys.argv[0])
    for smi in [
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]",
        "Cn1cc(NC=O)cc1C(=O)Nc1cc(C(=O)Nc2cc(C(=O)NCCC(N)=[NH2+])n(C)c2)n(C)c1",
        "OC(c1ccncc1)c1ccc(OCC[NH+]2CCCC2)cc1",
    ]:
        print(scorer.score(smi))
