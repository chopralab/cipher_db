from collections.abc import Mapping
from datetime import datetime
from multiprocessing import Process
import sys

import bson
import pymongo

from pymongo.collection import Collection

import json
import os
from multiprocessing import Process
from gridfs import GridFS
sys.path.insert(0, '../utils/')
from auto_mine_pubchem import PubChem_Miner
from get_scores import ScoreFetcher
from get_tree import TreeBuilder
from vis_tree import TreeVisualizer
import bson
from datetime import datetime
import urllib

sys.path.insert(0, '/Users/ludo/src/CANDO-master')
import cando as cnd

# We will change the login system later
login = open('../utils/login.txt','r')
username = login.readline().replace('\n','')
username = urllib.parse.quote(username)
password = login.readline().replace('\n','')
password = urllib.parse.quote(password)

#-------------------------------
# ------ UPDATE HELPERS --------
#-------------------------------

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

def update_reactivity_synthscore(id, result, coll):
    filled = []
    empty = []
    update = []
    for key,val in result.items():
        if key == 'smiles':
            filled.append((key,val))
        elif key == 'synthscores':
            print(key)
            print(val.items())
            for prop_key,prop_val in val.items():
                if prop_val is None or prop_val == '':
                    empty.append(prop_key)

    print(filled)
    print(empty)

    scorer = ScoreFetcher()

    for elem in filled:
        smi = elem[1]
        for data in empty:
            if data == 'sascore':
                output = scorer.get_SAScore(smi)
                update.append((data,output))
            elif data == 'scscore':
                output = scorer.get_SCScore(smi)
                update.append((data,output))

    print(update)

    # Perform database update
    for elem in update:
        key = elem[0]
        val = elem[1]
        print(key)
        print(val)
        doc = coll.find_one_and_update(
            {"_id" : id},
            {"$set":
                {f"synthscores.{key}": val}
            },upsert=True
        )

    # Perform date and time update
    doc = coll.find_one_and_update(
        {"_id" : id},
        {"$currentDate": {'modified': True}},
        upsert=True
    )

def update_reactivity_tree_viz(id, coll, tree, smi, db, num):
    tv = TreeVisualizer()
    graph = tv.tree_to_file(tree, smi)
    with open(f'{smi}.png', 'rb') as f:
        contents = f.read()
    fs = GridFS(db)
    tree_id = fs.put(contents, filename=smi)
    doc = coll.find_one_and_update(
        {"_id" : id},
        {"$set":
            {f"askcos.images.{num}.imageID": tree_id,
             f"askcos.images.{num}.name": smi}
        },upsert=True
    )
    os.remove(f'{smi}.png')
    os.remove(smi)

def update_reactivity_tree_build(id, result, coll, db):
    filled = []
    update = []
    for key,val in result.items():
        if key == 'smiles':
            filled.append((key,val))
    
    if len(filled) != 1:
        print('Error processing SMILES data')
        return
    print(filled)

    smi = filled[0][1]

    # Replace with ASKCOS Host
    # NOTE(degraff): this should work now
    HOST = 'http://35.188.119.21'

    # Create tree builder instance and get parameters
    tb = TreeBuilder(HOST)
    with open('tree_params.json') as f:
        params = json.load(f)
    
    # Build the tree(s) and append the results
    js = tb.build_tree(smi, params)
    for key,val in js.items():
        update.append((key, val))

    # Perform database update
    for elem in update:
        key = elem[0]
        val = elem[1]
        print(key)
        print(val)
        doc = coll.find_one_and_update(
            {"_id" : id},
            {"$set":
                {f"askcos.{key}": val}
            },upsert=True
        )

    if 'trees' in js:
        counter = 0
        for tree in js['trees']:
            path = f'tree{counter}:{smi}'.replace('/', '%2F')
            update_reactivity_tree_viz(id, coll, tree, path, db, str(counter))
            counter += 1

    # Perform date and time update
    doc = coll.find_one_and_update(
        {"_id" : id},
        {"$currentDate": {'modified': True}},
        upsert=True
    )

def update_reactivity(id, db, coll):
    try:
        result = coll.find_one(id)
    except Exception as e:
        result = None
        print("Look up Failed!")
        print(e)
    update_reactivity_synthscore(id, result, coll)
    update_reactivity_tree_build(id, result, coll, db)


def update_binding_bandock(id, coll):
    mongo_login = 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
    client = pymongo.MongoClient(mongo_login)
    general = client['cipher_aspire']['general']
    try:
        orig = general.find_one(id)
    except Exception as e:
        result = None
        print("Look up Failed!")
        print(e)
        return
    try:
        result = coll.find_one(id)
    except:
        print("Document does not exists in this collection.")
        coll.insert(id)

    smi = orig['smiles'] 
    name = orig['name']
    doc = coll.find_one_and_update(
        {'_id': id},
        {"$set":
            {f'name': f'{name}', f'smiles': f'{smi}'}
        },upsert=True
    )
    sig = cnd.generate_signature_smi(smi, fp="rd_ecfp4", vect="int", dist="dice", 
            org="aspire", bs="coach", c_cutoff=0.0, p_cutoff=0.0, 
            percentile_cutoff=0.0, i_score="dCxP", save_sig=False)
    for i in sig.index:
        doc = coll.find_one_and_update(
            {'_id': id},
            {"$set":
                {f'BANDOCK.{i}': f'{sig.loc[i,0]}'}
            },upsert=True
        )

#-------------------------
# ------ TRIGGERS --------
#-------------------------

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
        print("Error Occurred")

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
        print("Error Occurred")

def reactivity_trigger():
    mongo_login = 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
    client = pymongo.MongoClient(mongo_login)
    reactivity = client['cipher_aspire']['reactivity']
    try:
        with reactivity.watch([{'$match': {'operationType': 'insert'}}]) as stream:
            for insert_change in stream:
                # Do something
                print(insert_change)
                print('-----')
                id = insert_change['documentKey']['_id']
                update_reactivity(id, client['cipher_aspire'], reactivity)
    except pymongo.errors.PyMongoError as e:
        # The ChangeStream encountered an unrecoverable error or the
        # resume attempt failed to recreate the cursor.
        print("Error Occurred")
        print(e)

def binding_bandock_trigger():
    mongo_login = 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
    client = pymongo.MongoClient(mongo_login)
    binding = client['cipher_aspire']['binding']
    general = client['cipher_aspire']['general']

    try:
        with general.watch([{'$match': {'operationType': 'update'}}]) as stream:
            for update_change in stream:
                # Do something
                print(update_change)
                print('-----')
                id = update_change['documentKey']['_id']
                update_binding_bandock(id, binding)
    except pymongo.errors.PyMongoError as e:
        # The ChangeStream encountered an unrecoverable error or the
        # resume attempt failed to recreate the cursor.
        print("Error Occurred")
        print(e)

def main():
    p1 = Process(target=general_trigger)
    p2 = Process(target=property_trigger)
    p3 = Process(target=reactivity_trigger)
    #p4 = Process(target=binding_bandock_trigger)
    p1.start()
    p2.start()
    p3.start()
    #p4.start()

if __name__ == '__main__':
    main()
