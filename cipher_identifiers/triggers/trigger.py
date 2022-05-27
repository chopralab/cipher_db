import os
import mongoengine as me
import pymongo
import sys
import argparse
from rdkit import Chem

sys.path.append("../../")
from cipher_identifiers.utils.compounds import id_compound_from_smiles

parser = argparse.ArgumentParser()
parser.add_argument("--testing", action="store_true")
args = parser.parse_args()

if args.testing:
    URI = os.environ["TESTING_URI"]
else:
    URI = os.environ["MONGO_URI"]

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
                inchikey = Chem.MolToInchiKey(m)
                if id != Chem.MolToInchiKey(m):
                    compounds_coll.delete_one({"_id": id})
                try:
                    id_compound_from_smiles(smiles, inchikey)
                except Exception as e:
                    print(e)
                    compounds_coll.delete_one({"_id": id})
    except Exception as e:
        print(e)


if __name__ == "__main__":
    compounds_trigger()
