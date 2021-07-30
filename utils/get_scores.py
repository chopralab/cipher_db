from standalone_model_numpy import SCScorer
# import RAscore_NN
import rdkit.Chem as Chem
import rdkit.Chem.Draw as Draw
import sascorer
import os
project_root = os.path.dirname(os.path.dirname(__file__))

class ScoreFetcher():
    def __init__(self):
        # Import the SCScorer numpy model
        self.SCScorer = SCScorer()
        self.SCScorer.restore(os.path.join(project_root, 'models', 'full_reaxys_model_2048bool', 'model.ckpt-10654.as_numpy.json.gz'), FP_len=2048)
        # self.RAScorer = RAscore_NN.RAScorerNN()

    def get_SAScore(self, smi: str):
        mol = Chem.MolFromSmiles(smi)
        sco = sascorer.calculateScore(mol)
        # print('SAScore: %.4f <--- %s' % (sco, smi))
        return sco

    def get_SCScore(self, smi: str):
        (smi, sco) = self.SCScorer.get_score_from_smi(smi)
        # print('SCScore: %.4f <--- %s' % (sco, smi))
        return sco

    # def get_RAScore(self, smi: str):
    #     sco = self.RAScorer.predict(smi)
    #     print('RAScore: %.4f <--- %s' % (sco, smi))
    #     return sco

if __name__ == '__main__':
    model = ScoreFetcher()
    # List of SMILES strings
    smis = ['CC1=C(COC2=CC(OC)=C(CN3CCCC[C@H]3C(O)=O)C(OC)=C2)C=CC=C1C4=CC=CC=C4',
            'COC1=CC(OCC2=CC=CC(=C2C)C2=CC=C3OCCOC3=C2)=CC(OC)=C1CN[C@H](CO)C(O)=O',
            'COC1=CC(OCC2=CC=CC(=C2C)C2=CC=CC=C2)=CC(OC)=C1CN[C@H](CC(F)(F)F)C1=CC=CC=C1',
            'O=C1CCC2C(C[C@@H](C)[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4(C#C)O)=C1',
            'COC1=CC(OCC2=C(C)C(=CC=C2)C2=CC=CC=C2)=CC(OC)=C1CN1CCC[C@@H](C1)C(O)=O',
            'CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5']
    for smi in smis:
        model.get_SCScore(smi)
        model.get_SAScore(smi)
        # model.get_RAScore(smi)
        print("\n")