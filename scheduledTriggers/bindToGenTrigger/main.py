import pymongo
import certifi
import requests
import datetime


def bindToGenTrigger(event, context):
    username = 'Fong'
    password = 'password12345'
    mongo_login = 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
    client = pymongo.MongoClient(mongo_login, tlsCAFile=certifi.where())
    db = client['cipher_aspire']
    binding = client['cipher_aspire']['binding']
    general = client['cipher_aspire']['general']


    for mol_data in db.binding.find():
        old_data = {
            "inchikey": mol_data.get("inchikey")
        }
        old_data = general.find_one(old_data)

        # update general collection based on its existence
        # update general collection based on datetime? (requires datetime stamp in binding collection)
        if not (old_data):
            inchiKey = mol_data.get("inchikey")
            name = mol_data.get("name")

            if (inchiKey):
                try:
                    request_string = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{}/{}/property/IUPACName,InChI,InChIKey,CanonicalSMILES/JSON".format("InChiKey", inchiKey)
                    response = requests.get(request_string)
                    response = response.json()["PropertyTable"]["Properties"][0]

                    temp = {
                        "name": mol_data.get("name"),
                        "IUPACname": response.get("IUPACName"),
                        "actualInChi": response.get("InChI"),
                        "inchikey": response.get("InChIKey"),
                        "smiles": response.get("CanonicalSMILES"),
                        "CID": response.get("CID"),
                        "date_modified": datetime.datetime.now().strftime("%m/%d/%Y, %H:%M:%S"),
                    }
                    general.insert_one(temp)
                    continue

                except:
                    print("WARNING: Unable to update using SMILES: ", mol_data.get("name"))

            if (name):
                request_string = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{}/{}/property/IUPACName,InChI,InChIKey,CanonicalSMILES/JSON".format("name", name)
                response = requests.get(request_string)
                try:
                    response = response.json()["PropertyTable"]["Properties"][0]

                    temp = {
                        "name": mol_data.get("name"),
                        "IUPACname": response.get("IUPACName"),
                        "actualInChi": response.get("InChI"),
                        "inchikey": response.get("InChIKey"),
                        "CID": response.get("CID"),
                        "smiles": response.get("CanonicalSMILES"),
                        "date_modified": datetime.datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
                    }
                    print(temp)
                    general.insert_one(temp)
                    continue

                except:
                    print("WARNING: Unable to update using name: ", mol_data.get("name"))
            else:
                print("WARNING: No smiles or name field in binding data")
                
