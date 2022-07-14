"""adapated from https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/sascorer.py"""
import math
import sys
from typing import List

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


class SAScorer:
    def __init__(
        self, fpscoress: List[List[float]], lower_bound: float = -4, upper_bound: float = 2.5
    ):
        self.fp_score = {fp: score for score, *fps in fpscoress for fp in fps}
        self.lb = lower_bound
        self.range = upper_bound - lower_bound

    def __call__(self, smi: str) -> float:
        return self.score(smi)

    def score(self, smi: str) -> float:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"invalid SMILES string! got: {smi}")

        fp = rdMolDescriptors.GetMorganFingerprint(mol, 2)
        fps = fp.GetNonzeroElements()

        score1 = 0.0
        nf = 0

        for bid, v in fps.items():
            nf += v
            score1 += v * self.fp_score.get(bid, -4)
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

        sascore = 11.0 - (sascore - self.lb + 1) / (self.range) * 9.0

        if sascore > 8:
            sascore = 8 + math.log(sascore - 8)

        return np.clip(sascore, 1, 10)


if __name__ == "__main__":
    scorer = SAScorer(sys.argv[0])
    for smi in [
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]",
        "Cn1cc(NC=O)cc1C(=O)Nc1cc(C(=O)Nc2cc(C(=O)NCCC(N)=[NH2+])n(C)c2)n(C)c1",
        "OC(c1ccncc1)c1ccc(OCC[NH+]2CCCC2)cc1",
    ]:
        print(scorer.score(smi))
