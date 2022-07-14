import sys
sys.path.append("../..")
from cipher_cando.triggers.trigger import (
    update_cando
)
from cipher_identifiers.docs.docs import (
    Compounds
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

    for i in Compounds.objects():
        inchikey = i.inchikey
        smi = i.smiles
        update_cando(inchikey, smi, "dCxP", "aspire")
