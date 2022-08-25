import sys
import cando as cnd

sys.path.append("../..")
from cipher_cando.triggers.trigger import (
    update_cando
)
from cipher_cando.docs.docs import (
    Cando
)

from cipher_identifiers.docs.docs import (
    Compounds,
    Biomolecules,
    Models
)

if __name__ == "__main__":
    import argparse
    import mongoengine as me
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("--testing", action="store_true", help="Whether or not to run on the test database")
    args = parser.parse_args()

    LOGIN = open("../../utils/login.txt", "r")
    USERNAME = LOGIN.readline().replace("\n", "")
    PASSWORD = LOGIN.readline().replace("\n", "")
    LOGIN.close()

    if args.testing:
        URI = (
            "mongodb+srv://"
            + USERNAME
            + ":"
            + PASSWORD
            + "@aspirecluster0.hmj3q.mongodb.net/cipher_testing?retryWrites=true&w=majority"
        )
    else:
        URI = (
            "mongodb+srv://"
            + USERNAME
            + ":"
            + PASSWORD
            + "@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority"
        )

    me.connect(host=URI)

    i_score = "dCxP"
    mid = Models.objects(__raw__={"source": 'CANDO', "parameters": { "version": f'{cnd.__version__}', "interaction_score": f'{i_score}' } } ).get().cipher_mid
    '''
    prot_ids = ["P14416","P35372","P35462","P41143","P41145","P41146","P42262","Q13224"]
    bm_mid = Models.objects(__raw__={"source": 'AlphaFold', f"parameters": {"version": 'v2.0'} }).get().cipher_mid
    biomols = []
    for p in prot_ids:
        biomols.append(Biomolecules.objects(__raw__={ "uniprot_id": p, "cipher_mid": bm_mid }).get().cipher_bmid)
    '''
    for i in Compounds.objects():
        inchikey = i.inchikey
        smi = i.smiles
        for b in biomols:
            if not Cando.objects(inchikey=inchikey, cipher_bmid=b, cipher_mid=mid):
                update_cando(inchikey, smi, i_score, "aspire")
                update_cando(inchikey, smi, i_score, "homo_sapien")
                break

