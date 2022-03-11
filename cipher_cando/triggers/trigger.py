import mongoengine as me
import pymongo
import sys
import argparse

sys.path.insert(0, '/Users/ludo/src/CANDO-master')
import cando as cnd

sys.path.append("../../")
from cipher_identifiers.utils.compounds import id_compound_from_smiles
from cipher_identifiers.docs.docs import (
    Models,
    gen_unique_model_id,
    Biomolecules,
    gen_unique_biomol_id,
    Binding_sites,
    gen_unique_bindingsite_id,
    Ligands
)

#from ciper_cando.utils.utils import #####
from cipher_cando.docs.docs import (
    Cando,
    Biosignatures,
    knn_tuple,
    KNN,
    gen_unique_cando_id,
    gen_unique_biosig_id
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
    id2name_map = {"P14416":"DRD2","P35372":"OPRM1","P35462":"DRD3","P41143":"OPRD1","P41145":"OPRK1","P41146":"OPRL1","P42262": "GRIA2","Q13224":"GRIN2B"}
    sig = cnd.generate_signature_smi(smi, fp="rd_ecfp4", vect="int", dist="dice",
            org=organism, bs="coach", c_cutoff=0.0, p_cutoff=0.0,
            percentile_cutoff=0.0, i_score=i_score, save_sig=False)
    sig_bs = cnd.generate_signature_smi(smi, fp="rd_ecfp4", vect="int", dist="dice",
            org=organism, bs="coach", c_cutoff=0.0, p_cutoff=0.0,
            percentile_cutoff=0.0, i_score=i_score, save_sig=False,
            lig_name=True)
        
    # Biomolecule Model id
    # TODO: Replace this hardcoded query with metadata from CANDO
    if not Models.objects(__raw__={"source": 'AlphaFold', "parameters": { "version": 'v2' } }):
        bm_mid = gen_unique_model_id()
        model = Models()
        model.cipher_mid = bm_mid
        model.source = 'AlphaFold'
        model.parameters = { "version" : 2 }
        model.save()
    else:
        bm_mid = Models.objects(__raw__={"source": 'AlphaFold', f"parameters": { "version" : 'v2'} }).get().cipher_mid
    # TODO: run function to get predicted binding sites used in CANDO for each interaction
    # Will populate the binding site collection and add info to the binding data below


    # Create a dictionary of biomolecule internal ids
    # insert new documents if biomolecule does not already
    # exist in Biomolecules collection
    bmids = {}
    for i in sig.index:
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
                biomolecule.name = id2name_map[i]
                biomolecule.uniprot_id = i
                biomolecule.pdb_id = ''
                biomolecule.chain_id = ''
                biomolecule.cipher_mid = bm_mid
            else:
                biomolecule.cipher_bmid = bmid
                biomolecule.name = id2name_map[i]
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
