import sys
import pymongo
from multiprocessing import Process
sys.path.insert(0, '../utils/')
from auto_mine_pubchem import PubChem_Miner
from get_scores import ScoreFetcher
import bson
from datetime import datetime

# We will change the login system later
login = open('login.txt','r')
username = login.readline().replace('\n','')
password = login.readline().replace('\n','')

#-------------------------------
# ------ UPDATE HELPERS --------
#-------------------------------

def update_general(id,coll):
    try:
        result = coll.find_one(id)
        filled = []
        empty = []
        update = []
        for key,val in result.items():
            if isinstance(val,bson.objectid.ObjectId):
                pass
            elif val is not None and val != '':
                filled.append((key,val))
            else:
                empty.append(key)

        print(filled)
        print(empty)

        for elem in filled:
            input_type = elem[0]
            input = elem[1]
            for data in empty:
                if data == 'name':
                    output = PubChem_Miner.get_compound_info(input_type,input,'synonyms',data,'TXT')
                    output = output.split('\n')[0]
                    update.append((data,output))
                else:
                    output = PubChem_Miner.get_compound_info(input_type,input,'property',data,'TXT')
                    update.append((data,output))

        print(update)

        # Perform mined PubChem database update
        for elem in update:
            key = elem[0]
            val = elem[1]
            doc = coll.find_one_and_update(
                {"_id" : id},
                {"$set":
                    {key: val}
                },upsert=True
            )

        # Perform data and time update
        doc = coll.find_one_and_update(
            {"_id" : id},
            {"$currentDate": {'modified': True}},
            upsert=True
        )
    except:
        result = None
        print("Look Up Failed!")

def update_property_pubchem(id,result,coll):
    filled = []
    empty = []
    update = []
    for key,val in result.items():
        if key == 'inchikey':
            filled.append((key,val))
        elif key == 'pubchem':
            print(key)
            print(val.items())
            for prop_key,prop_val in val.items():
                if prop_val is not None and prop_val != '':
                    pass
                else:
                    empty.append(prop_key)

    print(filled)
    print(empty)

    for elem in filled:
        input_type = elem[0]
        input = elem[1]
        for data in empty:
            if data == 'name':
               output = PubChem_Miner.get_compound_info(input_type,input,'synonyms',data,'TXT')
               output = output.split('\n')[0]
               update.append((data,output))
            else:
                output = PubChem_Miner.get_compound_info(input_type,input,'property',data,'TXT')
                update.append((data,output))

    print(update)
    # Perform mined PubChem database update
    for elem in update:
        key = elem[0]
        val = elem[1]
        doc = coll.find_one_and_update(
            {"_id" : id, "database.dbname" : "pubchem"},
            {"$set":
                {key: val}
            },upsert=True
        )
    
    # Perform data and time update
    doc = coll.find_one_and_update(
        {"_id" : id, "database.dbname" : "pubchem"},
        {"$currentDate": {'modified': True}},
        upsert=True
    )

def update_property(id, coll):
    try:
        result = coll.find_one(id)
    except:
        result = None
        print("Look Up Failed!")
    
    update_property_pubchem(id,result,coll)

def update_reactivity_sascore(id, result, coll):
    result = coll.find_one(id)
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
                if prop_val is not None and prop_val != '':
                    pass
                else:
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

def update_reactivity(id, coll):
    try:
        result = coll.find_one(id)
    except Exception as e:
        result = None
        print("Look up Failed!")
        print(e)
    
    update_reactivity_sascore(id, result, coll)

#-------------------------
# ------ TRIGGERS --------
#-------------------------

def general_trigger():
    mongo_login = 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
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
    mongo_login = 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
    client = pymongo.MongoClient(mongo_login)
    properties = client['cipher_aspire']['properties']
    try:
        with properties.watch([{'$match': {'operationType': 'insert'}}]) as stream:
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
                update_reactivity(id, reactivity)
    except pymongo.errors.PyMongoError as e:
        # The ChangeStream encountered an unrecoverable error or the
        # resume attempt failed to recreate the cursor.
        print("Error Occurred")
        print(e)

def main():
    p1 = Process(target=general_trigger)
    p2 = Process(target=property_trigger)
    p3 = Process(target=reactivity_trigger)
    p1.start()
    p2.start()
    p3.start()

if __name__ == '__main__':
    main()