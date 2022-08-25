import sys
sys.path.append("../..")
import math
from sklearn.metrics import pairwise_distances
from scipy import stats

from cipher_identifiers.docs.docs import (
    Compounds,
    Biomolecules,
    Models
)
from cipher_cando.docs.docs import (
    Cando,
    Biosignatures,
    knn_tuple,
    Knn,
    gen_unique_knn_id
)

def update_knn(biosig_id, k):
    biosig = Biosignatures.objects(cipher_sig_id=biosig_id).get()
    all_scores = []
    all_cmpds = []
    all_smis = []
    for cmpd in Compounds.objects():
        all_cmpds.append(cmpd.inchikey)
        all_smis.append(cmpd.smiles)
        temp_scores = []
        for receptor in biosig.receptors:
            temp_scores.append(Cando.objects(inchikey=cmpd.inchikey, cipher_bmid=receptor,cipher_mid=biosig.cipher_mid).get().interaction_score)
        all_scores.append(temp_scores)
    distances = pairwise_distances([biosig.scores], all_scores, metric='cosine')[0]
    all_ranks = stats.rankdata(distances, method='max')
    knns = list(zip(all_cmpds, all_smis, distances, all_ranks))
    knns = sorted(knns, key=lambda x: x[3] if not math.isnan(x[2]) else 100000)

    if not Knn.objects(cipher_sig_id=biosig.cipher_sig_id):
        kid = gen_unique_knn_id()
    else:
        kid = Knn.objects(cipher_sig_id=biosig.cipher_sig_id).get().cipher_knn_id

    knn_tuples = []
    for i in knns[:k]:
        knn_tuples.append(knn_tuple(inchikey = i[0], smiles = i[1], cosine_dist = i[2], rank = i[3]))
    
    knn = Knn()
    knn.cipher_knn_id = kid
    knn.cipher_sig_id = biosig.cipher_sig_id
    knn.neighbors = knn_tuples
    knn.save()


if __name__ == "__main__":
    import argparse
    import mongoengine as me
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("--testing", action="store_true", help="Whether or not to run on the test database")
    parser.add_argument("--biosig", type=str, help="cipher_bsid corresponding to the biosignature to which all compounds will be compared")
    parser.add_argument("--k", type=str, help="Top k ranked compounds when compared to the biosignature, Default: 10")
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
    if args.k:
        k = int(args.k)
    else:
        k = 10
    if args.biosig:
        update_knn(args.biosig,k)
    else:
        # Update KNN with PZM21 as desired biosig
        # Hardcoded for now
        update_knn("AN2GZO",k)

