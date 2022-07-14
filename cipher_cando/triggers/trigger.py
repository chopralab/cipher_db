import mongoengine as me
import pymongo
import sys
import argparse
import math

import cando as cnd
from sklearn.metrics import pairwise_distances
from scipy import stats

sys.path.append("../../")
from cipher_identifiers.utils.compounds import id_compound_from_smiles
from cipher_identifiers.docs.docs import (
    Models,
    gen_unique_model_id,
    Biomolecules,
    gen_unique_biomol_id,
    Binding_sites,
    gen_unique_bindingsite_id,
    Ligands,
    Compounds
)

#from ciper_cando.utils.utils import #####
from cipher_cando.docs.docs import (
    Cando,
    Biosignatures,
    knn_tuple,
    Knn,
    gen_unique_cando_id,
    gen_unique_biosig_id,
    gen_unique_knn_id
)

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

def update_biosig(biomols, i_score, scores, description=''):
    if len(biomols) != len(scores):
        raise me.ValidationError("Length of desired dynamical biosignature does not match the number of receptors in biosignature.")
    mid = Models.objects(__raw__={"source": 'CANDO', "parameters": { "version": f'{cnd.__version__}', "interaction_score": f'{i_score}' } }).get().cipher_mid
    if not Biosignatures.objects(cipher_mid=mid, receptors__all=biomols):
        biosig_id = gen_unique_biosig_id()
    else:
        biosig_id = Biosignatures.objects(cipher_mid=mid, receptors__all=biomols).get().cipher_sig_id
    biosig = Biosignatures()
    biosig.cipher_sig_id = biosig_id
    biosig.cipher_mid = mid
    biosig.receptors = biomols
    biosig.scores = scores
    biosig.description = description
    biosig.save()

    return biosig.cipher_sig_id

def update_knn(biosig_id):
    biosig = Biosignatures.objects(cipher_sig_id=biosig_id).get()
    all_scores = []
    all_cmpds = []
    all_smis = []
    for cmpd in Compounds.objects()[:10]:
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
    for k in knns:
        knn_tuples.append(knn_tuple(inchikey = k[0], smiles = k[1], cosine_dist = k[2], rank = k[3]))
    
    knn = Knn()
    knn.cipher_knn_id = kid
    knn.cipher_sig_id = biosig.cipher_sig_id
    knn.neighbors = knn_tuples
    knn.save()


def update_cando(inchikey, smi, i_score, organism):
    # CANDO model id
    if not Models.objects(__raw__={"source": 'CANDO', "parameters": { "version": f'{cnd.__version__}', "interaction_score": f'{i_score}' } } ):
        mid = gen_unique_model_id()
        model = Models()
        model.cipher_mid = mid
        model.source = 'CANDO'
        model.parameters = { "version" : f'{cnd.__version__}', 'interaction_score' : f'{i_score}'}
        model.save()
    else:
        mid = Models.objects(__raw__={"source": 'CANDO', f"parameters": { "version" : f'{cnd.__version__}', 'interaction_score' : f'{i_score}'}}).get().cipher_mid
      
    # TODO: Replace with dataframe mapping with CANDO
    #id2name_map = {"P14416":"DRD2","P35372":"OPRM1","P35462":"DRD3","P41143":"OPRD1","P41145":"OPRK1","P41146":"OPRL1","P42262":"GRIA2","Q13224":"GRIN2B"}
    #id2name_map = {"P14416":"DRD2","P14416-F2":"D2SDR","P35372":"OPRM1","P35462":"DRD3","P41143":"OPRD1","P41145":"OPRK1","P41146":"OPRL1","P42262":"GRIA2","Q13224":"GRIN2B"}
    sig = cnd.generate_signature_smi(smi, fp="rd_ecfp4", vect="int", dist="dice",
            org=organism, bs="coach", c_cutoff=0.0, p_cutoff=0.0,
            percentile_cutoff=0.0, i_score=i_score, save_sig=False)
    sig_bs = cnd.generate_signature_smi(smi, fp="rd_ecfp4", vect="int", dist="dice",
            org=organism, bs="coach", c_cutoff=0.0, p_cutoff=0.0,
            percentile_cutoff=0.0, i_score=i_score, save_sig=False,
            lig_name=True)
        
    # Create a dictionary of biomolecule internal ids
    # insert new documents if biomolecule does not already
    # exist in Biomolecules collection
    bmids = {}
    for i in sig.index:
        bm_info = cnd.get_prot_info(i,org=organism)
        
        # Biomolecule Model id
        if not Models.objects(__raw__={"source": bm_info['method'], "parameters": {"version": bm_info['method_version']} }):
            bm_mid = gen_unique_model_id()
            model = Models()
            model.cipher_mid = bm_mid
            model.source = bm_info['method']
            model.parameters = { "version" : bm_info['method_version'] }
            model.save()
        else:
            bm_mid = Models.objects(__raw__={"source": bm_info['method'], f"parameters": {"version": bm_info['method_version']} }).get().cipher_mid

        if Biomolecules.objects(__raw__={ "uniprot_id": i, "cipher_mid": bm_mid }):
            bmids[i] = Biomolecules.objects(__raw__={ "uniprot_id": i, "cipher_mid": bm_mid }).get().cipher_bmid
        elif Biomolecules.objects(__raw__={ "pdb_id": f'{i[:4]}', "chain_id": f'{i[4:]}', "cipher_mid": bm_mid }):
            bmids[i] = Biomolecules.objects(__raw__={ "pdb_id": f'{i[:4]}', "chain_id": f'{i[4:]}', "cipher_mid": bm_mid }).get().cipher_bmid
        else:
            # Use CANDO to get all metadata for biomolecule and insert the new document
            # using the hardcoded dictionary for now
            bmid = gen_unique_biomol_id()
            biomolecule = Biomolecules()
            if len(i) > 5:
                biomolecule.cipher_bmid = bmid
                biomolecule.name = bm_info['uniprotRecommendedName']
                #biomolecule.name = id2name_map[i]
                biomolecule.uniprot_id = i
                biomolecule.pdb_id = ''
                biomolecule.chain_id = ''
                biomolecule.cipher_mid = bm_mid
            else:
                biomolecule.cipher_bmid = bmid
                biomolecule.name = bm_info['uniprotRecommendedName']
                biomolecule.uniprot_id = ''
                biomolecule.pdb_id = i[:4]
                biomolecule.chain_id = i[4:]
                biomolecule.cipher_mid = bm_mid
            biomolecule.save()
            bmids[i] = bmid

    # Binding site for compound-protein interaction
    bsids = {}
    for i in sig.index:
        if not Binding_sites.objects(__raw__={"cipher_bmid": bmids[i], "template_id": sig_bs.loc[i,0]}):
            bsid = gen_unique_bindingsite_id()
            binding_site = Binding_sites()
            binding_site.cipher_bsid = bsid
            binding_site.cipher_bmid = bmids[i]
            binding_site.template_id = sig_bs.loc[i,0]
            binding_site.save()
            bsids[i] = bsid
        else:
            bsids[i] = Binding_sites.objects(__raw__={"cipher_bmid": bmids[i], "template_id": sig_bs.loc[i,0]}).get().cipher_bsid

    for i in sig.index:
        # Push new document
        if not Cando.objects(__raw__={"inchikey": inchikey, "cipher_mid": mid, "cipher_bmid": bmids[i]}):
            cando = Cando()
            cando.cipher_cando_id = gen_unique_cando_id() 
            cando.inchikey = inchikey
            cando.smiles = smi
            cando.cipher_mid = mid
            cando.cipher_bmid = bmids[i]
            cando.cipher_bsid = bsids[i]
            cando.interaction_score = sig.loc[i,0]
            cando.save()

def cando_trigger():
    if args.testing:
        compounds_coll = MONGO_CLIENT.cipher_testing.compounds
    else:
        compounds_coll = MONGO_CLIENT.cipher_aspire.compounds
    try:
        with compounds_coll.watch([{"$match": {"operationType": "insert"}}]) as stream:
            for change in stream:
                d = change["fullDocument"]
                inchikey = d["inchikey"]
                smi = d["smiles"]

                update_cando(inchikey, smi, "dCxP", "aspire")
    except Exception as e:
        print(e)


if __name__ == "__main__":
    cando_trigger()
