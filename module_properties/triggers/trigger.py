import mongoengine as me
import pymongo
import sys
import argparse

sys.path.append("../../")
from module_properties.utils.properties import insert_properties_from_smiles

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


def properties_trigger():
    if args.testing:
        compounds_coll = MONGO_CLIENT.cipher_testing.compounds
    else:
        compounds_coll = MONGO_CLIENT.cipher_aspire.compounds

    try:
        with compounds_coll.watch([{"$match": {"operationType": "insert"}}]) as stream:
            for change in stream:
                doc = change["fullDocument"]
                smiles = doc["smiles"]
                try:
                    insert_properties_from_smiles(smiles)
                except Exception as e:
                    print(e)
    except Exception as e:
        print(e)


if __name__ == "__main__":
    properties_trigger()
