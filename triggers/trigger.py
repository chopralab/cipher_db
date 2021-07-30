from collections.abc import Mapping
from datetime import datetime
from multiprocessing import Process
import sys

import bson
import pymongo
from pymongo.collection import Collection

sys.path.insert(0, 'utils/')
from auto_mine_pubchem import PubChem_Miner

# We will change the login system later
with open('login.txt','r') as login:
    username = login.readline().replace('\n','')
    password = login.readline().replace('\n','')

def update_general(id: Mapping, collection: Collection):
    try:
        result = collection.find_one(id)

        filled = []
        empty = []
        for key, val in result.items():
            if isinstance(val, bson.objectid.ObjectId):
                pass
            elif val: # is not None and val != '':
                filled.append((key, val))
            else:
                empty.append(key)

        print(filled)
        print(empty)

        update = []
        for input_type, input in filled:
            for data in empty:
                if data == 'name':
                    output = PubChem_Miner.get_compound_info(
                        input_type, input, 'synonyms', data, 'TXT'
                    )
                    output = output.split('\n')[0]
                    update.append((data, output))
                else:
                    output = PubChem_Miner.get_compound_info(
                        input_type, input, 'property', data, 'TXT'
                    )
                    update.append((data, output))

        print(update)

        # Perform mined PubChem database update
        for key, val in update:
            doc = collection.find_one_and_update(
                {"_id" : id},
                {"$set": {key: val}},
                upsert=True
            )

        # Perform data and time update
        doc = collection.find_one_and_update(
            {"_id" : id},
            {"$currentDate": {'modified': True}},
            upsert=True
        )
    except:
        result = None
        print("Look Up Failed!")

def update_property_pubchem(
        id, result: Mapping, collection: Collection
    ):
    """update the database entry in the collection with specified id with result

    Parameters
    ----------
    id : [type]
        the id of the entry to update in the collection
    result : Mapping
        the entry corresponding to id
    collection : Collection
        the collection in which to upate the entry
    """
    filled = []
    data_to_update = []  # empty (previously)
    for key, val in result.items():
        if key == 'inchikey':
            filled.append((key, val))
        elif key == 'pubchem':
            print(val.items())
            for prop_key, prop_val in val.items():
                if prop_val: # is not None and prop_val != '':
                    pass
                else:
                    data_to_update.append(prop_key)

    print(filled)
    print(data_to_update)

    update = []
    for input_type, input in filled:
        # input_type is always 'inchikey'
        for data in data_to_update:
            if data == 'name':
                output = PubChem_Miner.get_compound_info(
                    input_type, input, 'synonyms', data, 'TXT'
                )
                output = output.split('\n')[0]
                update.append((data, output))
            else:
                output = PubChem_Miner.get_compound_info(
                    input_type, input, 'property', data, 'TXT'
                )
                update.append((data, output))

    print(update)
    for key, val in update:
        doc = collection.find_one_and_update(
            {"_id" : id, "database.dbname" : "pubchem"},
            {"$set": {key: val}},
            upsert=True
        )
    
    doc = collection.find_one_and_update(
        {"_id" : id, "database.dbname" : "pubchem"},
        {"$currentDate": {'modified': True}},
        upsert=True
    )

def update_property(id, collection: Collection):
    try:
        result = collection.find_one(id)
    except:
        result = None
        print("Look Up Failed!")
    
    update_property_pubchem(id, result, collection)

def general_trigger():
    mongo_login = f'mongodb+srv://{username}:{password}@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
    client = pymongo.MongoClient(mongo_login)
    db = client['cipher_aspire']
    general = client['cipher_aspire']['general']
    try:
        with general.watch([{'$match': {'operationType': 'insert'}}]) as stream:
            for insert_change in stream:
                # Do something
                print(insert_change)
                id = insert_change['documentKey']['_id']
                update_general(id, general)
    except pymongo.errors.PyMongoError:
        # The ChangeStream encountered an unrecoverable error or the
        # resume attempt failed to recreate the cursor.
        print("Error Occured")

def property_trigger():
    mongo_login = f'mongodb+srv://{username}:{password}@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
    client = pymongo.MongoClient(mongo_login)

    properties = client['cipher_aspire']['properties']
    try:
        with (
            properties.watch([{'$match': {'operationType': 'insert'}}])
        ) as stream:
            for insert_change in stream:
                # Do something
                print(insert_change)
                id = insert_change['documentKey']['_id']
                update_property(id, properties)
    except pymongo.errors.PyMongoError:
        # The ChangeStream encountered an unrecoverable error or the
        # resume attempt failed to recreate the cursor.
        print("Error Occured")

def main():
    p1 = Process(target=general_trigger)
    p2 = Process(target=property_trigger)
    p1.start()
    p2.start()

if __name__ == '__main__':
    main()