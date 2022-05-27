import pymongo as pmg
import mongoengine as me
import os
import argparse
import sys
sys.path.append("../")
from cipher_identifiers.utils.compounds import id_compound_from_smiles


parser = argparse.ArgumentParser()
parser.add_argument("--testing", action="store_true", help="Flag if inserting into the testing database")
args = parser.parse_args()

if args.testing:
    URI = os.environ["TESTING_URI"]
else:
    URI = os.environ["MONGO_URI"]

me.connect(host=URI)
MONGO_CLIENT = pmg.MongoClient(os.environ["MONGO_URI"])

if args.testing:
    compounds_coll = MONGO_CLIENT.cipher_testing.compounds
else:
    compounds_coll = MONGO_CLIENT.cipher_aspire.compounds

cursor = compounds_coll.find({})
for doc in cursor:
    id = doc["_id"]
    remove_op = compounds_coll.delete_one({"_id": id})
    print(remove_op)


