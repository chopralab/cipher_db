import json
import requests
import datetime

if __name__ == "__main__":
    import sys
    sys.path.append("../../")

from rdkit import Chem
from rdkit.Chem.Descriptors import (
    ExactMolWt,
    HeavyAtomMolWt,
    MolWt,
    MaxPartialCharge,
    MinPartialCharge,
    NumRadicalElectrons,
    NumValenceElectrons,
)
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.QED import default, weights_max, weights_none
from rdkit.Chem.Lipinski import (
    NHOHCount,
    NOCount,
    NumHAcceptors,
    NumHDonors,
    NumHeteroatoms,
    NumRotatableBonds,
    RingCount,
)
from cipher_properties.utils.errors import (
    InvalidJsonError,
    StructNotFoundError,
    InChiKeyNotFoundError,
    InvalidPubChemPropertyError,
    InvalidRequestError,
    InvalidSMILESError,
)
from cipher_identifiers.docs.docs import Compounds
from cipher_properties.docs.docs import Properties, Pubchem, RDKit


def __format_request_url(inchikey, selected_properties):
    """
    Formats the PugREST request URL for mining chemical property information from the PubChem database

    Parameters
    ----------
    properties: list, required
        The list which contains the specified chemical property information

    Returns
    -------
    url: string
        The URL to call to mine the chemical property information from PubChem

    Raises
    ------
    InvalidPubChemPropertyError
        If one or more of the chemical properties provided are invalid based on the PubChem RESTful API property table
        See https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest for more details
    """

    full_property_list = [
        "MolecularFormula",
        "MolecularWeight",
        "XLogP",
        "ExactMass",
        "MonoisotopicMass",
        "TPSA",
        "Complexity",
        "Charge",
        "HBondDonorCount",
        "HBondAcceptorCount",
        "RotatableBondCount",
        "HeavyAtomCount",
        "IsotopeAtomCount",
        "AtomStereoCount",
        "DefinedAtomStereoCount",
        "UndefinedAtomStereoCount",
        "BondStereoCount",
        "DefinedBondStereoCount",
        "UndefinedBondStereoCount",
        "CovalentUnitCount",
        "Volume3D",
        "XStericQuadrupole3D",
        "YStericQuadrupole3D",
        "ZStericQuadrupole3D",
        "FeatureCount3D",
        "FeatureAcceptorCount3D",
        "FeatureDonorCount3D",
        "FeatureAnionCount3D",
        "FeatureCationCount3D",
        "FeatureRingCount3D",
        "FeatureHydrophobeCount3D",
        "ConformerModelRMSD3D",
        "EffectiveRotorCount3D",
        "ConformerCount3D",
        "Fingerprint2D",
    ]

    identifier_list = [
        "CanonicalSMILES",
        "IsomericSMILES",
        "InChI",
        "InChIKey",
        "IUPACName",
        "Title",
    ]

    url_base = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/" + inchikey + "/property/"
    )
    url_middle = ""
    url_postfix = "/json"

    if len(selected_properties) == 0:
        raise InvalidPubChemPropertyError(
            "No properties selected - Please select at least one property"
        )
    elif "All" in selected_properties:
        selected_properties = full_property_list
        for prop in selected_properties:
            url_middle += prop
            url_middle += ","
    else:
        for prop in selected_properties:
            if prop in identifier_list:
                raise InvalidPubChemPropertyError(
                    "Identifier selected - selected properties cannot include CanonicalSMILES, IsomericSMILES, InChI, InChIKey, IUPACName, Title"
                )
            elif prop not in full_property_list:
                raise InvalidPubChemPropertyError(
                    "Invalid property selected - Please check the PubChem RESTful API property table to see a list of valid properites"
                )
            else:
                url_middle += prop
                url_middle += ","

    url_middle = url_middle[:-1]

    url = url_base + url_middle + url_postfix
    return url


def __make_url_request(inchikey, url):
    """
    Gets the associated JSON data from the formatted PubChem URL request

    Parameters
    ----------
    inchikey: string, required
        The InChI Key of the compound to make the URL request for

    Raises
    ------
    InvalidRequestError
        If the request does not return a status code in the 200's (i.e. status code relating to an error)
    """
    response = requests.get(url)
    if response:
        data = json.loads(response.text)
        return data
    else:
        return None


def insert_properties_from_smiles(smiles, inchikey, properties=["All"]):
    try:
        m = Chem.MolFromSmiles(smiles)
    except:
        raise InvalidSMILESError(
            "Invalid SMILES - The provided SMILES "
            + smiles
            + " cannot be converted into a valid molecule"
        )

    prop = Properties()
    prop.inchikey = inchikey
    prop.smiles = smiles
    prop.validate()
    
    try:
        url = __format_request_url(inchikey, properties)
        data = __make_url_request(inchikey, url)
        if data is not None:
            data = data["PropertyTable"]["Properties"][0]
        else:
            data = {}

        pc = Pubchem()
        for key, val in data.items():
            setattr(pc, key, val)
        pc.modified = datetime.datetime.utcnow
    except Exception as e:
        print(e)

    rd = RDKit()
    rd.cipher_mid = "002"
    rd.MolWt = MolWt(m)
    rd.ExactMolWt = ExactMolWt(m)
    rd.HeavyAtomMolWt = HeavyAtomMolWt(m)
    rd.MaxPartialCharge = MaxPartialCharge(m)
    rd.MinPartialCharge = MinPartialCharge(m)
    rd.NumRadicalElectrons = NumRadicalElectrons(m)
    rd.NumValenceElectrons = NumValenceElectrons(m)
    rd.MolLogP = MolLogP(m)
    rd.MaxQED = weights_max(m)
    rd.MeanQED = default(m)
    rd.NoneQED = weights_none(m)
    rd.NHOHCount = NHOHCount(m)
    rd.NOCount = NOCount(m)
    rd.NumHAcceptors = NumHAcceptors(m)
    rd.NumHDonors = NumHDonors(m)
    rd.NumHeteroatoms = NumHeteroatoms(m)
    rd.NumRotateableBonds = NumRotatableBonds(m)
    rd.RingCount = RingCount(m)
    rd.modified = datetime.datetime.utcnow

    prop.pubchem = pc
    prop.rdkit = rd
    prop.save()
    print("Inserted properties for object with InChI Key {} to the Properties Collectino of the database".format(inchikey))


if __name__ == "__main__":
    import argparse
    import os
    import mongoengine as me

    parser = argparse.ArgumentParser()
    parser.add_argument("--testing", action="store_true", help="Flag if inserting into the testing database")
    parser.add_argument("--update", action="store_true")
    parser.add_argument("--count", action="store_true")
    args = parser.parse_args()

    if args.testing:
        URI = os.environ["TESTING_URI"]
    else:
        URI = os.environ["MONGO_URI"]

    me.connect(host=URI)

    if args.update:
        for comp in Compounds.objects:
            try:
                smiles = comp.smiles
                inchikey = comp.inchikey
                if Properties.objects.with_id(inchikey) is None:
                     insert_properties_from_smiles(smiles, inchikey)
            except Exception as e:
                print(e)
    if args.count:
        print(Properties.objects().count())
    else:
        pass