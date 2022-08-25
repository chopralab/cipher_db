import mongoengine as me
import pymongo
import sys
import argparse
import math, time

import cando as cnd
from sklearn.metrics import pairwise_distances
from scipy import stats
import pandas as pd

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
        if i[0].isalpha() and pd.isna(bm_info['PDB_ID']):
            bm_info['GENE_NAME']=" "
            bm_info['METHOD']=" "
        elif i[0].isnumeric() and pd.isna(bm_info['UNIPROT_ID']):
            bm_info['GENE_NAME']=" "
            bm_info['UNIPROT_ID']=" "
            bm_info['METHOD']="solved"
        # Biomolecule Model id
        if not Models.objects(__raw__={"source": bm_info['METHOD'], "parameters": {"version": bm_info['METHOD_VERSION']} }):
            bm_mid = gen_unique_model_id()
            model = Models()
            model.cipher_mid = bm_mid
            model.source = bm_info['METHOD']
            model.parameters = { "version" : bm_info['METHOD_VERSION'] }
            model.save()
        else:
            bm_mid = Models.objects(__raw__={"source": bm_info['METHOD'], f"parameters": {"version": bm_info['METHOD_VERSION']} }).get().cipher_mid

        if Biomolecules.objects(__raw__={ "uniprot_id": i, "cipher_mid": bm_mid }):
            bmids[i] = Biomolecules.objects(__raw__={ "uniprot_id": i, "cipher_mid": bm_mid }).get().cipher_bmid
        elif Biomolecules.objects(__raw__={ "pdb_id": f'{i[:4]}', "chain_id": f'{i[4:]}', "cipher_mid": bm_mid }):
            bmids[i] = Biomolecules.objects(__raw__={ "pdb_id": f'{i[:4]}', "chain_id": f'{i[4:]}', "cipher_mid": bm_mid }).get().cipher_bmid
        else:
            # Use CANDO to get all metadata for biomolecule and insert the new document
            # TODO: Add other metadata fields to document and pull from metadata 
            bmid = gen_unique_biomol_id()
            biomolecule = Biomolecules()
            biomolecule.cipher_bmid = bmid
            biomolecule.name = bm_info['GENE_NAME']
            biomolecule.cipher_mid = bm_mid
            if i[0].isalpha():
                biomolecule.uniprot_id = i
                biomolecule.pdb_id = ''
                biomolecule.chain_id = ''
            else:
                biomolecule.uniprot_id = bm_info['UNIPROT_ID']
                biomolecule.pdb_id = i[:4]
                biomolecule.chain_id = i[4:]
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
                time.sleep(3)
                d = change["fullDocument"]
                inchikey = d["_id"]
                print(inchikey)
                smi = d["smiles"]
                print(smi)

                update_cando(inchikey, smi, "dCxP", "aspire")
                update_cando(inchikey, smi, "dCxP", "homo_sapien")
    except Exception as e:
        print(e)


if __name__ == "__main__":
    cando_trigger()
