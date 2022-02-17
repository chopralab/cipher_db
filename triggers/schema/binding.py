from gc import collect
from xml.dom.minidom import Element
import pymongo
import json
import requests
import pandas as pd
import deepdiff


class InvalidJsonError(Exception):
    pass


class InvalidFileError(Exception):
    pass

class InChiKeyNotFoundError(Exception):
    pass


class StructNotFoundError(Exception):
    pass


class pubchem():
    '''
    Class for interacting with PubChem binding information in the binding collection of the database

    Methods
    -------
    def insert_from_file()
    '''

    @staticmethod
    def __request_id_from_cid(cid):
        request_string = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{}/{}/property/CanonicalSMILES,InchiKey/JSON".format("cid", cid)
        response = requests.get(request_string)
        if response:
            response = json.loads(response.text)
            inchikey = response["PropertyTable"]["Properties"][0]["InChIKey"]
            return inchikey
        else:
            print("Request ID failed using CID, Status Code: " + str(response.status_code))
            return None

    @staticmethod
    def __make_hashable(unhashable_list):
        for i in range(len(unhashable_list)):
            for j in range(len(unhashable_list[i])):
                if type(unhashable_list[i][j]) == type([]):
                    new_tup = list(unhashable_list[i])
                    new_tup[j] = tuple(new_tup[j])
                    unhashable_list[i] = tuple(new_tup)
        return unhashable_list

    @staticmethod
    def insert_from_json_file(filename, source, receptor, collection):
        with open(filename) as json_file:
            data = json.load(json_file)
            for elem in data:
                inchikey = pubchem.__request_id_from_cid(elem['cid'])
                if collection.find_one({"inchikey": inchikey}) is not None:
                    if collection.find_one({"inchikey": inchikey, "pubchem."+receptor: {"$exists": True}}):
                        same_source = False
                        for doc in collection.find({"inchikey": inchikey}):
                            for key in doc["pubchem"][receptor]:
                                if key.startswith(source):
                                    same_source = True
                        if same_source:
                            identical = False
                            docs = collection.find({"inchikey": inchikey, "pubchem."+receptor+"."+source: {"$exists": True}})
                            for doc in docs:
                                if elem == doc["pubchem"][receptor][source]:
                                    identical = True
                                else:
                                    elem_set = set(pubchem.__make_hashable(list(elem.items())))
                                    doc_set = set(pubchem.__make_hashable(list(doc["pubchem"][receptor][source].items())))
                                    diff = dict(elem_set - doc_set)
                                    if 'units' not in diff.keys():
                                        identical = True
                                    else:
                                        if source in list(doc["pubchem"][receptor].keys()):
                                            units = doc["pubchem"][receptor][source]["units"]
                                            collection.find_one_and_update({"inchikey": inchikey}, {"$rename": {"pubchem."+receptor+"."+source: "pubchem."+receptor+"."+source+"_"+units}})
                            if identical:
                                print("Identical entry already exists in the database for " + inchikey)
                            else:
                                units = elem["units"]
                                collection.find_one_and_update(
                                    {"inchikey": inchikey},
                                    {"$set": {"pubchem."+receptor +"."+source+"_"+units: elem}}
                                )
                                collection.find_one_and_update(
                                    {"inchikey": inchikey},
                                    {"$currentDate": {'modified': True}},
                                    upsert=True
                                )
                        else:
                            collection.find_one_and_update(
                                {"inchikey": inchikey},
                                {"$set": {"pubchem."+receptor+"."+source: elem}}
                            )
                            collection.find_one_and_update(
                                {"inchikey": inchikey},
                                {"$currentDate": {'modified': True}},
                                upsert=True
                            )
                    else:
                        collection.find_one_and_update({"inchikey": inchikey}, {"$set": {"pubchem."+receptor+"."+source: elem}})
                else:
                    entry = {"inchikey": inchikey, "pubchem": {receptor: {source: elem}}}
                    db_entry = collection.insert_one(entry)
                    collection.find_one_and_update(
                        {"inchikey": inchikey},
                        {"$currentDate": {'modified': True}},
                        upsert=True
                    )
                    print(db_entry.inserted_id)

    @staticmethod
    def insert(inchikey):
        '''
        '''
        pass

    @staticmethod
    def remove(inchikey):
        '''
        '''
        pass

    @staticmethod
    def edit(inchikey, json):
        '''
        '''
        pass


class cando():
    '''
    Class for interacting with CANDO binding data in the binding collection of the database

    Methods
    -------
    insert(inchikey)
        Runs the CANDO module on the specified compound and updates the binding collection of the database with the resulting information
    remove(inchikey)
        Removes the CANDO struct of the specified compound from the binding collection of the database
    edit(inchikey, json)
        Edits the CANDO struct of the specified compound with the information in the specified JSON file 
    '''
    @staticmethod
    def insert(inchikey):
        '''
        '''
        pass

    @staticmethod
    def remove(inchikey):
        '''
        '''
        pass

    @staticmethod
    def edit(inchikey, json):
        '''
        '''
        pass


class desi_exp():
    '''
    Class for interacting with DESI experiments in the binding collection of the database

    Methods
    -------
    insert_from_tsv_file(fname, collection, cipher_eid)
        Method for inserting compounds from the following tab seperated values (tsv) file into the database
    '''

    @staticmethod
    def insert_from_tsv_file(fname, collection, cipher_eid):
        '''
        Method for inserting compounds from the following tab seperated values (tsv) file into the database

        Parameters
        ----------
        fname: string, required
            The name of the tsv file to be parsed
        collection: Mongo DB collection, required
            The database collection to insert the information into (unless testing use binding)
        cipher_eid: string, required
            The the cipher experiment id that corresponds to the specific desi experiment which genearted the datafile
        '''
        df = pd.read_csv(fname, sep="\t")
        for index, row in df.iterrows():
            inchikey = row["inchikey"]
            receptor = row["receptor"]
            data = {
                    "cipher_eid": cipher_eid,
                    "units": row["units"],
                    "compound_concentration": row["compound_concentration"],
                    "response": row["response"],
                    "receptor_concentration": row["receptor_concentration"],
                    "competitor": row["competitor"],
                    "competitor_concentration": row["competitor_concentration"],
                    "control": row["control"],
                    "ratio": row["ratio"]
            }
            if collection.find_one({"inchikey": inchikey}) is not None:
                if collection.find_one({"inchikey": inchikey, "desi."+receptor: {"$exists": True}}):
                    experiments = collection.find_one({"inchikey": inchikey})["desi"][receptor]
                    duplicate = False
                    for e in experiments:
                        if e == data:
                            duplicate = True
                    if duplicate:
                        print("Duplicate Assay: " + str(data))
                    else:
                        collection.find_one_and_update(
                            {"inchikey": inchikey},
                            {"$push":
                                    {"desi."+receptor: data}
                            }
                        )
                        collection.find_one_and_update(
                            {"inchikey": inchikey},
                            {"$currentDate": {'desi.modified': True}} ,
                            upsert=True
                        )
                else:
                    collection.find_one_and_update(
                        {"inchikey": inchikey},
                        {"$set": {
                            "desi."+receptor:[data]
                            }
                        }
                    )
            else:
                entry = {
                    "inchikey": inchikey, 
                    "desi": {
                        receptor: [data]
                        }
                    }
                db_entry = collection.insert_one(entry)
                collection.find_one_and_update(
                   {"inchikey": inchikey},
                   {"$currentDate": {'desi.modified': True}} ,
                   upsert=True
                )
                print(db_entry.inserted_id)

    @staticmethod
    def insert(json):
        '''
        '''
        pass

    @staticmethod
    def remove(inchikey):
        '''
        '''
        pass

    @staticmethod
    def edit(inchikey, json):
        '''
        '''
        pass
