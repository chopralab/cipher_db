import sys
sys.path.append("../..")
from cipher_cando.triggers.trigger import (
    update_knn,
    update_biosig
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

    # Define Biosig --> PZM21
    # i.e., receptors and scores
    # Hard code list of proteins, but extract corresponding bmid
    prot_ids = ["P14416","P35372","P35462","P41143","P41145","P41146","P42262","Q13224"]
    bm_mid = Models.objects(__raw__={"source": 'AlphaFold', f"parameters": {"version": 'v2.0'} }).get().cipher_mid
    biomols = []
    for p in prot_ids:
        if Biomolecules.objects(__raw__={ "uniprot_id": p, "cipher_mid": bm_mid }):
            biomols.append(Biomolecules.objects(__raw__={ "uniprot_id": p, "cipher_mid": bm_mid }).get().cipher_bmid)
    # hard code scores...for now
    scores = [0.660,0.560,0.670,0.570,0.580,0.570,0.979,0.379]
    biosig_id = update_biosig(biomols, "dCxP", scores, description="PZM21 signature using ASPIRE protein library")
    
    # Update KNN with PZM21 as desired biosig
    update_knn(biosig_id)

