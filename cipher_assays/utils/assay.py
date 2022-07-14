import json

from numpy import insert, source
import shortuuid
import requests
from time import sleep

if __name__ == "__main__":
    import sys
    sys.path.append("../../")

from cipher_assays.docs.docs import Assays
from cipher_assays.docs.docs import gen_unique_assay_id

from cipher_identifiers.docs.docs import Compounds
from cipher_identifiers.utils.compounds import id_compound_from_smiles, id_compounds_from_inchikey

def __get_info_from_cid(cid):
    url ="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/property/InChIKey,CanonicalSMILES/JSON".format(cid)
    response = requests.get(url)
    if response:
        smiles = json.loads(response.content)["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
        inchikey = json.loads(response.content)["PropertyTable"]["Properties"][0]["InChIKey"]
        sleep(0.25)
        return smiles, inchikey
    else:
        raise requests.exceptions.HTTPError("Status code {}".format(response.status_code))

def insert_pubchem_assays_from_json(fname, receptor, source):
    with open(fname, 'r') as f:
        data = json.load(f)
        for doc in data:
            try:
                cid = doc["cid"]
            except:
                print("CID field not found in documnet")
                continue
            try:
                smiles, inchikey = __get_info_from_cid(cid)
                if Compounds.objects.with_id(inchikey) is None:
                    id_compounds_from_inchikey(inchikey)
                assay = Assays()
                assay.cipher_aid = gen_unique_assay_id()
                assay.inchikey = inchikey
                assay.smiles = smiles
                assay.source = source
                assay.receptor = receptor
                for key, val in doc.items():
                    setattr(assay, key, val)
                assay.save()
                print("Inserted Assay with Cipher AID {} into the Assays Collection of the Database".format(assay.cipher_aid))
            except Exception as e:
                print(e)



def insert_bulk_assays():
     # RECEPTOR CONSTANTS
    DELTA = "deltaOR"
    KAPPA = "kappaOR"
    MU = "muOR"
    NOCICEPTIN = "nociceptinOR"
    SIGMA = "sigmaOR"
    DRD2 = "DRD2"
    DRD3 = "DRD3"
    AMPAR = "AMPAR"
    NMDAR = "NMDAR"

    # SOURCE CONSTANTS
    CHEMBL = "chembl"
    DRUGBANK = "drugbank"
    IUPHAR = "iuphar"
    
    # DELTA
    chembl_delta_or_1 = "../data/PUBCHEM/DELTA/delta_opioid_receptor_chembl_drugs_1.json"
    chembl_delta_or_2 = "../data/PUBCHEM/DELTA/delta_opioid_receptor_chembl_drugs_2.json"
    drugbank_delta_or = "../data/PUBCHEM/DELTA/delta_opioid_receptor_drugbank_drugs.json"
    iuphar_delta_or = "../data/PUBCHEM/DELTA/delta_opioid_receptor_pharmacological_ligands.json"

    # KAPPA
    chembl_kappa_or_1 = "../data/PUBCHEM/KAPPA/kappa_opioid_receptor_chembl_drugs.json"
    chembl_kappa_or_2 = "../data/PUBCHEM/KAPPA/kappa_opioid_receptor_chembl_drugs_2.json"
    drugbank_kappa_or = "../data/PUBCHEM/KAPPA/kappa_opioid_receptor_drugbank_drugs.json"
    iuphar_kappa_or = "../data/PUBCHEM/KAPPA/kappa_opioid_receptor_IUPHAR_ligands.json"

    # MU
    chembl_mu_or = "../data/PUBCHEM/MU/pubchem/mu_opioid_receptor_chembl_drugs.json"
    drugbank_mu_or = "../data/PUBCHEM/MU/pubchem/mu_opioid_receptor_drugbank_drugs.json"
    iuphar_mu_or = "../data/PUBCHEM/MU/pubchem/mu_opioid_receptor_iuphar_ligands.json"    

    # NOCICEPTIN
    chembl_nociceptin_or = "../data/PUBCHEM/NOCICEPTIN/nociceptin_opioid_receptor_chembl_drugs.json"
    drugbank_nociceptin_or = "../data/PUBCHEM/NOCICEPTIN/nociceptin_opioid_receptor_drugbank_drugs.json"
    iuphar_nociceptin_or = "../data/PUBCHEM/NOCICEPTIN/nociceptin_opioid_receptor_IUPHAR_ligands.json"

    # SIGMA
    chembl_sigma_or = "../data/PUBCHEM/SIGMA/sigma_non_opioid_receptor_chembl_drugs.json"
    drugbank_sigma_or = "../data/PUBCHEM/SIGMA/sigma_non_opioid_receptor_drugbank_drugs.json"
    iuphar_sigma_or = "../data/PUBCHEM/SIGMA/sigma_non_opioid_receptor_IUPHAR_ligands.json"
    
    # DRD2
    chembl_DRD2 = "../data/PUBCHEM/DRD2/DRD2_receptor_chembl_drugs.json"
    drugbank_DRD2 = "../data/PUBCHEM/DRD2/DRD2_receptor_drugbank_drugs.json"
    iuphar_DRD2 = "../data/PUBCHEM/DRD2/DRD2_receptor_iuphar_ligands.json"

    # DRD3
    chembl_DRD3 = "../data/PUBCHEM/DRD3/DRD3_receptor_chembl_drugs.json"
    drugbank_DRD3 = "../data/PUBCHEM/DRD3/DRD3_receptor_drugbank_drugs.json"
    iuphar_DRD3 = "../data/PUBCHEM/DRD3/DRD3_receptor_iuphar_ligands.json"

    # NDMAR
    chembl_NMDAR = "../data/PUBCHEM/NMDAR/NMDA_receptor_chembl_drugs.json"
    drugbank_NMDAR = "../data/PUBCHEM/NMDAR/NMDA_receptor_drugbank_drugs.json"
    iuphar_NMDAR = "../data/PUBCHEM/NMDAR/NMDA_receptor_iuphar_ligands.json"

    # AMPAR
    chembl_AMPAR = "../data/PUBCHEM/AMPAR/AMPA_receptor_chembl_drugs.json"
    drugbank_AMPAR = "../data/PUBCHEM/AMPAR/AMPA_receptor_drugbank_drugs.json"
    iuphar_AMPAR = "../data/PUBCHEM/AMPAR/AMPA_receptor_iuphar_ligands.json"

    print(exists(chembl_delta_or_1))
    print(exists(chembl_delta_or_2))
    print(exists(drugbank_delta_or))
    print(exists(iuphar_delta_or))

    print(exists(chembl_kappa_or_1))
    print(exists(chembl_kappa_or_2))
    print(exists(drugbank_kappa_or))
    print(exists(iuphar_kappa_or))

    print(exists(chembl_mu_or))
    print(exists(drugbank_mu_or))
    print(exists(iuphar_mu_or))

    print(exists(chembl_nociceptin_or))
    print(exists(drugbank_nociceptin_or))
    print(exists(iuphar_nociceptin_or))

    print(exists(chembl_sigma_or))
    print(exists(drugbank_sigma_or))
    print(exists(iuphar_sigma_or))

    print(exists(chembl_DRD2))
    print(exists(drugbank_DRD2))
    print(exists(iuphar_DRD2))

    print(exists(chembl_DRD3))
    print(exists(drugbank_DRD3))
    print(exists(iuphar_DRD3))

    print(exists(chembl_NMDAR))
    print(exists(drugbank_NMDAR))
    print(exists(iuphar_NMDAR))

    print(exists(chembl_AMPAR))
    print(exists(drugbank_AMPAR))
    print(exists(iuphar_AMPAR))

    
    insert_pubchem_assays_from_json(chembl_delta_or_1, receptor=DELTA, source=CHEMBL)
    insert_pubchem_assays_from_json(chembl_delta_or_2, receptor=DELTA, source=CHEMBL)
    insert_pubchem_assays_from_json(drugbank_delta_or, receptor=DELTA, source=DRUGBANK)
    insert_pubchem_assays_from_json(iuphar_delta_or, receptor=DELTA, source=IUPHAR)

    insert_pubchem_assays_from_json(chembl_kappa_or_1, receptor=KAPPA, source=CHEMBL)
    insert_pubchem_assays_from_json(chembl_kappa_or_2, receptor=KAPPA, source=CHEMBL)
    insert_pubchem_assays_from_json(drugbank_kappa_or, receptor=KAPPA, source=DRUGBANK)
    insert_pubchem_assays_from_json(iuphar_kappa_or, receptor=KAPPA, source=IUPHAR)

    insert_pubchem_assays_from_json(chembl_mu_or, receptor=MU, source=CHEMBL)
    insert_pubchem_assays_from_json(drugbank_mu_or, receptor=MU, source=DRUGBANK)
    insert_pubchem_assays_from_json(iuphar_mu_or, receptor=MU, source=IUPHAR)

    insert_pubchem_assays_from_json(chembl_nociceptin_or, receptor=NOCICEPTIN, source=CHEMBL)
    insert_pubchem_assays_from_json(drugbank_nociceptin_or, receptor=NOCICEPTIN, source=DRUGBANK)
    insert_pubchem_assays_from_json(iuphar_nociceptin_or, receptor=NOCICEPTIN, source=IUPHAR)

    insert_pubchem_assays_from_json(chembl_sigma_or, receptor=SIGMA, source=CHEMBL)
    insert_pubchem_assays_from_json(drugbank_sigma_or, receptor=SIGMA, source=DRUGBANK)
    insert_pubchem_assays_from_json(iuphar_sigma_or, receptor=SIGMA, source=IUPHAR)
    
    insert_pubchem_assays_from_json(chembl_DRD2, receptor=DRD2, source=CHEMBL)
    insert_pubchem_assays_from_json(drugbank_DRD2, receptor=DRD2, source=DRUGBANK)
    insert_pubchem_assays_from_json(iuphar_DRD2, receptor=DRD2, source=IUPHAR)

    insert_pubchem_assays_from_json(chembl_DRD3, receptor=DRD3, source=CHEMBL)
    insert_pubchem_assays_from_json(drugbank_DRD3, receptor=DRD3, source=DRUGBANK)
    insert_pubchem_assays_from_json(iuphar_DRD3, receptor=DRD3, source=IUPHAR)

    insert_pubchem_assays_from_json(chembl_NMDAR, receptor=NMDAR, source=CHEMBL)
    insert_pubchem_assays_from_json(drugbank_NMDAR, receptor=NMDAR, source=DRUGBANK)
    insert_pubchem_assays_from_json(iuphar_NMDAR, receptor=NMDAR, source=IUPHAR)

    insert_pubchem_assays_from_json(chembl_AMPAR, receptor=AMPAR, source=CHEMBL)
    insert_pubchem_assays_from_json(drugbank_AMPAR, receptor=AMPAR, source=DRUGBANK)
    insert_pubchem_assays_from_json(iuphar_AMPAR, receptor=AMPAR, source=IUPHAR)

if __name__ == "__main__":
    import argparse
    import os
    from os.path import exists
    import mongoengine as me

    # RECEPTOR CONSTANTS
    DELTA = "deltaOR"
    KAPPA = "kappaOR"
    MU = "muOR"
    NOCICEPTIN = "nociceptinOR"
    SIGMA = "sigmaOR"
    DRD2 = "DRD2"
    DRD3 = "DRD3"
    AMPAR = "AMPAR"
    NMDAR = "NMDAR"

    # SOURCE CONSTANTS
    CHEMBL = "chembl"
    DRUGBANK = "drugbank"
    IUPHAR = "iuphar"
    BIOASSAY = "pubchem bioassay"

    parser = argparse.ArgumentParser()
    parser.add_argument("--testing", action="store_true", help="Flag if inserting into the testing database")
    args = parser.parse_args()

    if args.testing:
        URI = os.environ["TESTING_URI"]
    else:
        URI = os.environ["MONGO_URI"]

    me.connect(host=URI)

    pzm21_bioassay_data = "../data/Bioassay/PZM21_BioAssay_Data.json"
    insert_pubchem_assays_from_json(pzm21_bioassay_data, receptor=MU, source=BIOASSAY)
