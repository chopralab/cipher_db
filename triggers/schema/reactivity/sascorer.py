"""adapated from https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/sascorer.py"""
import math
import pickle
import sys

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


class SAScorer:
    def __init__(self, fp_scores_pkl):
        with open(fp_scores_pkl, "rb") as fid:
            xss = pickle.load(fid)
            
        self.fscores = {x: float(xs[0]) for xs in xss for x in xs[1:]}

    def __call__(self, smi: str) -> float:
        return self.score(smi)

    def score(self, smi: str):
        mol = Chem.MolFromSmiles(smi)
        fp = rdMolDescriptors.GetMorganFingerprint(mol, 2)
        fps = fp.GetNonzeroElements()

        score1 = 0.0
        nf = 0
        for bid, v in fps.items():
            nf += v
            score1 += v * self.fscores.get(bid, -4)
        score1 /= nf

        num_atoms = mol.GetNumAtoms()
        num_chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        num_bridgeheads = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
        num_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        num_macrocycles = 0

        for x in mol.GetRingInfo().AtomRings():
            if len(x) > 8:
                num_macrocycles += 1

        sizePenalty = num_atoms**1.005 - num_atoms
        stereoPenalty = math.log10(num_chiral_centers + 1)
        bridgePenalty = math.log10(num_bridgeheads + 1)
        spiroPenalty = math.log10(num_spiro + 1)
        macrocyclePenalty = 0.0

        if num_macrocycles > 0:
            macrocyclePenalty = math.log10(2)

        score2 = -(sizePenalty + stereoPenalty + spiroPenalty + bridgePenalty + macrocyclePenalty)

        score3 = 0.0
        if num_atoms > len(fps):
            score3 = math.log(float(num_atoms) / len(fps)) * 0.5

        sascore = score1 + score2 + score3

        min = -4.0
        max = 2.5
        sascore = 11.0 - (sascore - min + 1) / (max - min) * 9.0
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
