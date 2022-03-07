from os import remove
import mongoengine as me
import pymongo
import sys
import argparse
from rdkit import Chem

sys.path.append("../../")
from module_identifiers.utils.compounds import id_compound_from_smiles

parser = argparse.ArgumentParser()
parser.add_argument("--testing", action="store_true")
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
MONGO_CLIENT = pymongo.MongoClient(URI)


def compounds_trigger():
    if args.testing:
        compounds_coll = MONGO_CLIENT.cipher_testing.compounds
    else:
        compounds_coll = MONGO_CLIENT.cipher_aspire.compounds
    try:
        with compounds_coll.watch([{"$match": {"operationType": "insert"}}]) as stream:
            for change in stream:
                doc = change["fullDocument"]
                smiles = doc["smiles"]
                id = change["documentKey"]["_id"]
                m = Chem.MolFromSmiles(smiles)
                if id != Chem.MolToInchiKey(m):
                    compounds_coll.delete_one({"_id": id})
                try:
                    id_compound_from_smiles(smiles)
                except Exception as e:
                    print(e)
                    compounds_coll.delete_one({"_id": id})
    except Exception as e:
        print(e)


if __name__ == "__main__":
    compounds_trigger()
