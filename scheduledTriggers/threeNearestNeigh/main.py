import torch
import numpy as np
import pymongo
import certifi

username = 'Fong'
password = 'password12345'
mongo_login = 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
client = pymongo.MongoClient(mongo_login, tlsCAFile=certifi.where())
db = client['cipher_aspire']

def threeNearestNeigh(event, context):
    cos = torch.nn.CosineSimilarity(dim=1, eps=1e-6)
    
    for mol_data in db.testing.find():

        # remove old dict
        three_nearest_dict_old = { "$unset": { "3NN" : mol_data.get("3NN") } }
        db.testing.update_one({"_id": mol_data.get("_id")}, three_nearest_dict_old)

        nearest_neigh_dict_old = mol_data.get("nearest_neigh")
        nearest_neigh_dict_old_unset = { "$unset": { "nearest_neigh" : mol_data.get("nearest_neigh") } }
        db.testing.update_one({"_id": mol_data.get("_id")}, nearest_neigh_dict_old_unset)

        # break
            
        nearest_neigh = []
        for mol_data_comp in db.testing.find():

            # fresh start for new molecule
            if ((mol_data.get("name") != mol_data_comp.get("name")) ):
                
                if (mol_data.get("CANDO")):

                    # get bandock score for mol
                    bandock_list = mol_data.get("CANDO").values()
                    score_list = torch.tensor(np.array([[float(info.get('BANDOCK_score')) for info in bandock_list]]))

                    # get bandock score for mol which are compared
                    bandock_list_comp = mol_data_comp.get("CANDO").values()
                    score_list_comp = torch.tensor(np.array([[float(info.get('BANDOCK_score')) for info in bandock_list_comp]]))

                    # calculate_similarity
                    output = cos(score_list, score_list_comp)
                    nearest_neigh.append({
                        "name": mol_data_comp.get("name"),
                        "similarity": output.cpu().detach().numpy()[0]
                    })
                    # print(nearest_neigh)

        # if not mol_data.get("3NN"):
        nearest_neigh = sorted(nearest_neigh, key = lambda i: i['similarity'],reverse=True) or nearest_neigh_dict_old
        three_nearest_dict = { "$set": { "3NN" : nearest_neigh[:3] } }
        db.testing.update_one({"_id": mol_data.get("_id")}, three_nearest_dict)

        nearest_neigh_dict_set = { "$set": { "nearest_neigh" : nearest_neigh } }
        db.testing.update_one({"_id": mol_data.get("_id")}, nearest_neigh_dict_set)
