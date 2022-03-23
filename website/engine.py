# Import statements
from unittest import result
from flask import Flask
from matplotlib import image
from flask_pymongo import PyMongo
from flask_mongoengine import MongoEngine
import mongoengine as me
import datetime
from gridfs import GridFS
import codecs
import base64
import pymongo
import json
import sys

sys.path.append("../")

from cipher_askcos.docs.docs import (
    Retrosynthesis,
    Difficulty,
    SyntheticTree,
    ChemicalNode,
    ReactionNode,
)
from cipher_assays.docs.docs import Assays
from cipher_cando.docs.docs import Cando, Biosignatures, KNN, knn_tuple
from cipher_identifiers.docs.docs import Compounds, Models, Biomolecules
from cipher_properties.docs.docs import Properties, RDKit, Pubchem

# ------------------------------------------------------
# --------------- Flask Initilization ------------------
# ------------------------------------------------------

# Uncomment this for logging
# logging.basicConfig(filename='logs/restful.log', level=logging.DEBUG)
# Set up Flask-Mongo DB connections
login = open("../utils/login.txt", "r")
username = login.readline().replace("\n", "")
password = login.readline().replace("\n", "")
login.close()
mongo_login = (
    "mongodb+srv://"
    + username
    + ":"
    + password
    + "@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority"
)

# Initilize PyMongo and MongoEngine
client = pymongo.MongoClient(mongo_login)
me.connect(host=mongo_login)


def return_compounds(identifier):
    """
    Returns chemical compound information for a compound specified by a specific identifier

    Parameters:
    -----------
    identifier: string, required
        identifer pertaining to the compound in question (Name, SMILES, InChI Key, Pubchem CID, InChI, IUPAC Name)

    Returns:
    --------
    A JSON with all identifying information on the compound loacted in the CIPHER Database or None if it is not found
    """

    if Compounds.objects(name=identifier).count() > 0:
        return json.loads(Compounds.objects(name=identifier).to_json())
    elif Compounds.objects(smiles=identifier).count() > 0:
        return json.loads(Compounds.objects(smiles=identifier).to_json())
    elif Compounds.objects(inchikey=identifier).count() > 0:
        return json.loads(Compounds.objects(inchikey=identifier).to_json())
    elif Compounds.objects(cid=identifier).count() > 0:
        return json.loads(Compounds.objects(cid=identifier).to_json())
    elif Compounds.objects(inchi=identifier).count() > 0:
        return json.loads(Compounds.objects(inchi=identifier).to_json())
    elif Compounds.objects(iupac=identifier).count() > 0:
        return json.loads(Compounds.objects(iupac=identifier).to_json())
    else:
        results = []
        for comp in Compounds.objects:
            if comp.synonyms is not None and isinstance(comp.synonyms, list):
                if identifier in comp.synonyms:
                    results.append(json.loads(comp.to_json()))
        return results


def return_properties(inchikey):
    """
    Returns chemical property information on the compound specified by the given InChI Key

    Parameters:
    -----------
    inchikey: string, required
        The InChI Key of the selected compound

    Returns:
    A JSON with the chemical property information for the selected compound store in the CIPHER database
    """
    prop = Properties.objects().with_id(inchikey)
    if prop is not None:
        return json.loads(prop.to_json())


def return_biosignature(inchikey):
    """
    Return the CANDO biosignature for the selected compound

    Parameters:
    -----------
    inchikey: string, required
        The InChi Key of the selected compound

    Returns:
    The CANDO biosignature of the selected compound in JSON format or None if not found
    """
    BIOSIG = [
        "muOR",
        "deltaOR",
        "kappaOR",
        "nociceptinOR",
        "DRD2",
        "D2LDR",
        "DRD3",
        "NMDAR",
        "AMPAR",
    ]
    BMIDS = {
        "1NN6YL": BIOSIG[0],
        "1W94ZT": BIOSIG[1],
        "2IMRXM": BIOSIG[2],
        "14LKRS": BIOSIG[3],
        "XXXXXX": BIOSIG[4],
        "24JXV0": BIOSIG[5],
        "1GKVIX": BIOSIG[6],
        "14SUK2": BIOSIG[7],
        "1OWFJU": BIOSIG[8],
    }
    CHEMBL = "chembl"
    DRUGBANK = "drugbank"
    IUPHAR = "iuphar"

    moa = {}
    binding = {}
    for assay in Assays.objects(inchikey=inchikey):
        if assay.source == CHEMBL:
            if assay.receptor not in moa:
                if hasattr(assay, "action"):
                    moa[assay.receptor] = assay.action.lower()
        if assay.source == DRUGBANK:
            if assay.receptor not in moa:
                if hasattr(assay, "drugaction"):
                    moa[assay.receptor] = assay.drugaction.lower()
        if assay.source == IUPHAR:
            if assay.receptor not in moa:
                if hasattr(assay, "action"):
                    moa[assay.receptor] = assay.action.lower()

    for receptor in BIOSIG:
        if receptor not in moa:
            moa[receptor] = "unknown effect"

    for entry in Cando.objects(inchikey=inchikey):
        if entry.cipher_bmid in BMIDS.keys():
            binding[BMIDS[entry.cipher_bmid]] = float(entry.interaction_score)

    for receptor in BMIDS.values():
        if receptor not in binding:
            binding[receptor] = -1.0

    return moa, binding


def return_assays(inchikey):
    """
    Returns binding assay information on the selected compound

    Parameters:
    -----------
    inchikey: string, required
        The InChI Key of the selected compound

    Returns:
    --------
    pubchem and desi binding assay information in JSON format or None if not found
    """
    if Assays.objects(inchikey=inchikey) is not None:
        return json.loads(Assays.objects(inchikey=inchikey).to_json())
    else:
        return []


def return_askcos_pathways(inchikey):
    """
    Returns all ASKCOS pathways that were found for the given molecule

    Parameters:
    -----------
    inchikey: string, required
        The InChi Key of the selected compound

    Returns:
    -------
    images: list(string)
        A list of base64 decoded image strings, see notes for how to convert the string to an image

    Notes:
    ------
    To decode the image use the following:
        decoded_image = open("path1.png", 'wb')
        decoded_image.write(base64.b64decode(images[0]))
        decoded_image.close()
    """
    retro = Retrosynthesis.objects().with_id(inchikey)
    images = []
    if retro is not None:
        if hasattr(retro, "trees"):
            trees = retro.trees
            for tree in trees:
                image = tree.image.read()
                image = codecs.encode(image, "base64").decode("utf-8")
                images.append(image)
    return images

def return_compound_image(inchikey):
    '''
    Returns the image of the compoud with the corresponding InChI Key

    Parameters:
    -----------
    inchikey: string, required
        The InChI Key of the compound to have its image retreived

    Returns:
    --------
    image: string
        The image of the corresponding compound in a utf-8 decoded string
    '''
    comp = Compounds.objects().with_id(inchikey)
    if comp is not None:
        if hasattr(comp, "image"):
            image = comp.image.read()
            image = codecs.encode(image, 'base64').decode('utf-8')
            return image
    return None

def return_desired_dynamic_biosignature():
    return return_biosignature("RMRJXGBAOAMLHD-IHFGGWKQSA-N")

def testing():
    compounds_id_info = return_compounds("RONZAEMNMFQXRA-UHFFFAOYSA-N")
    properties = return_properties("RONZAEMNMFQXRA-UHFFFAOYSA-N")
    moa, binding = return_biosignature("RONZAEMNMFQXRA-UHFFFAOYSA-N")
    assay = return_assays("RONZAEMNMFQXRA-UHFFFAOYSA-N")
    image = return_compound_image("RONZAEMNMFQXRA-UHFFFAOYSA-N")
    retro = return_askcos_pathways("CQOJHAJWCDJEAT-UHFFFAOYSA-N")
    print("Compound Identifying Information")
    print(compounds_id_info)
    print(" ")
    print("Compound Properties")
    print(properties)
    print(" ")
    print("MOA Vector")
    print(moa)
    print(" ")
    print("CANDO Binding Vector")
    print(binding)
    print(" ")
    print("Assay Data")
    print(assay)

    d_moa, d_binding = return_desired_dynamic_biosignature()

    print(" ")
    print("Desired Biosignature MOA Vector")
    print(d_moa)
    print(" ")
    print("Desired Biosignature CANDO Binding Vector")
    print(d_binding)

    decoded_image = open("temp/RONZAEMNMFQXRA-UHFFFAOYSA-N.svg", 'wb')
    decoded_image.write(base64.b64decode(image))
    decoded_image.close()

    for i in range(len(retro)):
        decoded_image = open("temp/RONZAEMNMFQXRA-UHFFFAOYSA-N_tree_{}.png".format(i), 'wb')
        decoded_image.write(base64.b64decode(retro[i]))
        decoded_image.close()


if __name__ == "__main__":
    testing()
