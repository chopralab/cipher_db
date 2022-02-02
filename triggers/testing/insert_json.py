# Import relevant packages
import pymongo
import json
import requests
import argparse

# Command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True)
parser.add_argument("--receptor", type=str, required=True)
parser.add_argument("--source", type=str, required=True)
args = parser.parse_args()

# Database login and client setup
username = "Fong"
password = "password12345"
mongo_login = 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
client = pymongo.MongoClient(mongo_login)
collection = client['cipher_aspire']['testing']

# Use PubChem restful API to request identifying information (inchikey, smiles)
def request_id_from_cid(cid):
    request_string = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{}/{}/property/CanonicalSMILES,InchiKey/JSON".format("cid", cid)
    response = requests.get(request_string)
    if response:
        response = json.loads(response.text)
        inchikey = response["PropertyTable"]["Properties"][0]["InChIKey"]
        return inchikey
    else:
        print("Request ID failed using CID, Status Code: " + str(response.status_code))
        return None

# Converts any list elements to tuples for a list of key-value pairs in a single nested JSON
def make_hashable(unhashable_list):
    for i in range(len(unhashable_list)):
        for j in range(len(unhashable_list[i])):
            if type(unhashable_list[i][j]) == type([]):
                new_tup = list(unhashable_list[i])
                new_tup[j] = tuple(new_tup[j])
                unhashable_list[i] = tuple(new_tup)
    return unhashable_list

# Add JSON information 
def import_binding_json():
    with open(args.filename) as json_file:
        data = json.load(json_file)
        for elem in data:
            smiles, inchikey = request_id_from_cid(elem['cid'])
            name = request_name_from_cid(elem['cid'])
            if collection.find_one({"inchikey":inchikey}) is not None:
                if collection.find_one({"inchikey":inchikey, args.receptor:{"$exists": True}}):
                    same_source = False
                    for doc in collection.find({"inchikey":inchikey}):
                            for key in doc[args.receptor]:
                                if key.startswith(args.source):
                                    same_source = True
                    if same_source:
                        identical = False
                        docs = collection.find({"inchikey":inchikey, "pubchem."+args.receptor+"."+args.source:{"$exists": True}})
                        for doc in docs:
                           if elem == doc[args.receptor][args.source]:
                               identical = True
                           else:
                               elem_set = set(make_hashable(list(elem.items())))
                               doc_set = set(make_hashable(list(doc[args.receptor][args.source].items())))
                               diff = dict(elem_set - doc_set)
                               if 'units' not in diff.keys():
                                   identical = True
                               else:
                                   if args.source in list(doc[args.receptor].keys()):
                                       units = doc[args.receptor][args.source]["units"]
                                       collection.find_one_and_update({"inchikey":inchikey}, {"$rename":{args.receptor+"."+args.source:args.receptor+"."+args.source+"_"+units}})
                        if identical:
                            print("Identical entry already exists in the database for " + inchikey)
                        else:
                            units = elem["units"]
                            collection.find_one_and_update(
                                {"inchikey":inchikey}, 
                                {"$set":{"pubchem."+args.receptor+"."+args.source+"_"+units:elem}}
                            )
                            collection.find_one_and_update(
                                {"inchikey":inchikey},
                                {"$currentDate": {'modified': True}},
                                upsert=True
                            )
                    else:
                        collection.find_one_and_update(
                            {"inchikey":inchikey}, 
                            {"$set":{"pubchem."+args.receptor+"."+args.source:elem}}
                        )
                        collection.find_one_and_update(
                            {"inchikey":inchikey},
                            {"$currentDate": {'modified': True}},
                            upsert=True
                        )
                else:
                    collection.find_one_and_update({"inchikey":inchikey}, {"$set":{"pubchem."+args.receptor+"."+args.source:elem}})
            else:
                entry = {"inchikey":inchikey,"pubchem":{args.receptor:{args.source:elem}}}
                db_entry = collection.insert_one(entry)
                collection.find_one_and_update(
                    {"inchikey":inchikey},
                    {"$currentDate": {'modified': True}},
                    upsert=True
                )
                print(db_entry.inserted_id)

import_binding_json()