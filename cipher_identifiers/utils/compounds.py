from ast import Str
from turtle import update
import requests
import json
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw
from time import sleep
import os
import codecs
import base64
RDLogger.DisableLog('rdApp.*')  

if __name__ == "__main__":
    import sys
    sys.path.append("../../")

# Import Custom Errors
from cipher_identifiers.utils.errors import (
    InvalidRequestError,
    CompoundNotFoundError,
    InvalidSMILESError,
)

# Import all validation functions
from cipher_identifiers.docs.docs import (
    validate_smiles,
    check_inchikey_in_compounds,
    check_mid_in_models,
    check_bmid_in_biomolecules,
    check_bsid_in_binding_sites,
)

# Import ME document
from cipher_identifiers.docs.docs import Compounds
from cipher_properties.docs.docs import Properties

def update_image(inchikey):
    comp = Compounds.objects().with_id(inchikey)
    if os.isdir("../temp/") is False:
        os.mkdir("../temp/")
    fname = "../temp/"+inchikey+".svg"
    mol = Chem.MolFromSmiles(comp.smiles)
    Draw.MolToFile(mol, fname)
    with open(fname, "rb") as f:
        comp.image.put(f, content_type="image/svg")
        comp.save()
    f.close()
    os.remove(fname)

def get_image_from_mol(mol, inchikey, comp):
    fname = "../temp/"+inchikey+".svg"
    Draw.MolToFile(mol, fname)
    with open(fname, "rb") as f:
        comp.image.put(f, content_type="image/svg")
    f.close()
    os.remove(fname)

def get_smiles_from_inchikey(inchikey):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{}/property/canonicalSMILES/JSON".format(inchikey)
    response = requests.get(url)
    if response:
        data = json.loads(response.text)
        smiles = data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
        return smiles
    else:
        raise InvalidRequestError(
            "Invalid Request - Request on InChi Key "
            + inchikey
            + " failed with status code "
            + str(response.status_code)
        )

def get_ids_from_inchikey(inchikey):
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
        + inchikey
        + "/property/IUPACName,Title/json"
    )
    response = requests.get(url)
    if response:
        data = json.loads(response.text)
        cid = data["PropertyTable"]["Properties"][0]["CID"]
        iupac = data["PropertyTable"]["Properties"][0]["IUPACName"]
        name = data["PropertyTable"]["Properties"][0]["Title"]
        return str(cid), iupac, name
    else:
        raise InvalidRequestError(
            "Invalid Request - Request on InChi Key "
            + inchikey
            + " failed with status code "
            + str(response.status_code)
        )

def get_synonyms_from_inchikey(inchikey):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{}/synonyms/JSON".format(inchikey)
    response = requests.get(url)
    if response:
        data = json.loads(response.text)
        synonyms = data["InformationList"]["Information"][0]["Synonym"]
        return synonyms
    else:
        raise InvalidRequestError(
            "Invalid Request - Request on CID "
            + inchikey
            + " failed with status code "
            + str(response.status_code)
        )

def id_compound_from_smiles(smiles, inchikey):
    comp = Compounds()
    try:
        m = Chem.MolFromSmiles(smiles)
    except:
        raise InvalidSMILESError(
            "Invalid SMILES - The provided SMILES "
            + smiles
            + " cannot be converted into a valid molecule"
        )

    comp.inchikey = inchikey
    comp.smiles = Chem.MolToSmiles(m, isomericSmiles=False)
    comp.inchi = Chem.MolToInchi(m)

    get_image_from_mol(m, inchikey, comp)

    try:
        comp.cid, comp.iupac, comp.name = get_ids_from_inchikey(inchikey)
        sleep(0.25)
        comp.synonyms = get_synonyms_from_inchikey(inchikey)
    except Exception as e:
        print(e)

    comp.save()
    print(
        "Inserted Object with InChI Key {} to Compounds Collection of the database".format(inchikey)
    )
    return comp.inchikey

def id_compounds_from_inchikey(inchikey):
    if Compounds.objects.with_id(inchikey) is not None:
        return False

    comp = Compounds()
    comp.inchikey = inchikey
    try:
        comp.smiles = get_smiles_from_inchikey(inchikey)
        sleep(0.25)
        comp.cid, comp.iupac, comp.name = get_ids_from_inchikey(inchikey)
        sleep(0.25)
        comp.synonyms = get_synonyms_from_inchikey(inchikey)
    except Exception as e:
        print(e)

    m = Chem.MolFromSmiles(comp.smiles)
    comp.inchi = Chem.MolToInchi(m)

    comp.save()
    print(
        "Inserted Object with InChI Key {} to Compounds Collection of the database".format(inchikey)
    )
    return True

def get_image_from_database(inchikey):
    comp = Compounds.objects().with_id(inchikey)
    if comp is not None:
        image = comp.image.read()
        image = codecs.encode(image, 'base64').decode('utf-8')
        decoded_image = open("../temp/{}.svg".format(inchikey), 'wb')
        decoded_image.write(base64.b64decode(image))
        decoded_image.close()


if __name__ == "__main__":
    import argparse
    import os
    import mongoengine as me

    parser = argparse.ArgumentParser()
    parser.add_argument("--testing", action="store_true", help="Flag if inserting into the testing database")
    parser.add_argument("--update", action="store_true", help="Flag for mass update opertion on the database")
    parser.add_argument("--smiles", type=str, help="The SMILES string to be added to the database")
    parser.add_argument("--inchikey", type=str, help="The SMILES string to be added to the database")
    parser.add_argument("--get-image", action="store_true")
    parser.add_argument("--count", action="store_true")
    parser.add_argument("--name", type=str)
    args = parser.parse_args()

    if args.testing:
        URI = os.environ["TESTING_URI"]
    else:
        URI = os.environ["MONGO_URI"]

    me.connect(host=URI)

    if args.update:
       for comp in Compounds.objects():
            try:
                update_image(comp.inchikey)
                print("Updated image for compound with InChI Key {}".format(comp.inchikey))
            except Exception as e:
                print(e)
    if args.get_image:
        if args.inchikey is not None:
            get_image_from_database(args.inchikey)
    if args.count:
        print(Compounds.objects().count())
    else:
        if args.smiles is not None and args.inchikey is not None:
            sleep(1)
            if args.name is not None:
                comp = Compounds(inchikey=args.inchikey,name=args.name,smiles=args.smiles)
                comp.save()