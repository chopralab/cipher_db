import imp
import pymongo
import json
import requests
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


class InvalidJsonError(Exception):
    pass


class StructNotFoundError(Exception):
    pass


class InChiKeyNotFoundError(Exception):
    pass


class InvalidPubChemPropertyError(Exception):
    pass


class InvalidRequestError(Exception):
    pass


class pubchem:
    """
    Class for interacting with the PubChem struct of the properties collection of the database

    Methods
    -------
    insert(inchikey)
        Mines physical/chemical property information on the specified compound from the PubChem database and inserts it into the properties collection of the database
    remove(inchikey)
        Removes the PubChem struct of the specified compound from the properties collection of the database
    """

    @staticmethod
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

    @staticmethod
    def __make_url_request(inchikey, url):
        """
        Gets the associated JSON data from the formatted PubChem URL request

        Parameters
        ----------
        inchikey: string, required
            The InChI Key of the compound to make the URL request for

        Returns
        -------
        data: dict
            The JSON formatted data from the PUG REST request

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
            raise InvalidRequestError(
                "Invalid Request - Request on InChi Key "
                + inchikey
                + " failed with status code "
                + response.status_code
            )

    @staticmethod
    def insert(inchikey, collection, compounds, properties=["All"]):
        """
        Mines physical/chemical property information on the specified compound from the PubChem database and inserts it into the properties collection of the database

        Parameters
        ----------
        inchikey: string, required
            The InChi Key of the compound to be mined
        collection: Mongo DB collection, required
            The Mongo DB collection where data operations will take place
        properties: list, default is ["All"]
            The list of chemical properties to mine from PubChem, set single list element to "All" to mine all possible properties
        Raises
        ------
        InChiKeyNotFoundError
            If the InChiKey does not exist in the compounds collection of the database
        """

        if compounds.find_one({"inchikey": inchikey}) is None:
            raise InChiKeyNotFoundError(
                "The provided compound's InChiKey was not found in the compounds colletion of the database"
            )

        url = pubchem.__format_request_url(inchikey, properties)
        data = pubchem.__make_url_request(inchikey, url)
        data = data["PropertyTable"]["Properties"][0]
        data.pop("CID")

        if collection.find_one({"inchikey": inchikey}) is not None:
            if collection.find_one({"inchikey": inchikey, "pubchem": {"$exists": True}}):
                prev_data = collection.find_one({"inchikey": inchikey})["pubchem"]
                prev_data.pop("modified")
                data = {**data, **prev_data}
            collection.find_one_and_update({"inchikey": inchikey}, {"$set": {"pubchem": data}})
            collection.find_one_and_update(
                {"inchikey": inchikey}, {"$currentDate": {"pubchem.modified": True}}
            )
        else:
            entry = {"inchikey": inchikey, "pubchem": data}
            db_entry = collection.insert_one(entry)
            collection.find_one_and_update(
                {"inchikey": inchikey}, {"$currentDate": {"pubchem.modified": True}}
            )
            print(db_entry.inserted_id)

    @staticmethod
    def remove(inchikey):
        """
        Removes the PubChem struct of the specified compound from the properties collection of the database

        Parameters
        ----------
        inchikey: string, required
            The InChi Key of the entry which has its PubChem struct deleted

        Raises
        ------
        InChiKeyNotFoundError
            If the InChiKey does not exist in the properties collection of the database
        StructNotFoundError
            If a PubChem struct is not found for the given InChi Key entry in the properties collection
        """
        pass


class rdkit:
    @staticmethod
    def __calc_properties(inchikey, compounds, cipher_mid):
        """
        Method for calculating properties using RDKit and fomatting the data entry for database insertion

        Parameters
        ----------
        inchikey: string, required
            The InChI Key of the compound which will have its properties calculated
        compounds: Mongo DB collection, required
            The compounds collection of the database for constraint enforcement
        cipher_mid: string, required
            The model ID of the version of RDKit used

        Returns
        -------
        data: dict
            The formatted data entry for RDKit filled with calculated properties and descriptors
        """
        smiles = compounds.find_one({"inchikey": inchikey})["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        data = {
            "cipher_mid": cipher_mid,
            "MolWt": MolWt(mol),
            "ExactMolWt": ExactMolWt(mol),
            "HeavyAtomMolWt": HeavyAtomMolWt(mol),
            "MaxPartialCharge": MaxPartialCharge(mol),
            "MinPartialCharge": MinPartialCharge(mol),
            "NumRadicalElectrons": NumRadicalElectrons(mol),
            "NumValenceElectrons": NumValenceElectrons(mol),
            "MolLogP": MolLogP(mol),
            "MaxQED": weights_max(mol),
            "MeanQED": default(mol),
            "NoneQED": weights_none(mol),
            "NHOHCount": NHOHCount(mol),
            "NOCount": NOCount(mol),
            "NumHAcceptors": NumHAcceptors(mol),
            "NumHDonors": NumHDonors(mol),
            "NumHeteroatoms": NumHeteroatoms(mol),
            "NumRotatableBonds": NumRotatableBonds(mol),
            "RingCount": RingCount(mol),
        }
        return data

    @staticmethod
    def insert(inchikey, collection, compounds, cipher_mid):
        """
        Method for inserting chemical properties and descriptors calculated by RDKit to the properties collection of the databases

        Parameters
        ----------
        inchikey: string, required
            The InChI of the compound which will have its properties calculated
        collection: Mongo DB collection, required
            The database collection to insert the information into (should be properties)
        compounds: Mongo DB collection, required
            The compounds collection of the database to check certian constraints, access SMILES string information
        cipher_mid: string, required
            The model ID in the models collection of the database that corresponds to the specific verion of RDKit used

        Raises
        ------
        InChiKeyNotFoundError
            If the InChIKey was not found in the compounds collection of the database or the SMILES information was not present
        """
        if compounds.find_one({"inchikey": inchikey}) is None:
            raise InChiKeyNotFoundError(
                "The provided compound's InChiKey was not found in the compounds colletion of the database"
            )

        if "smiles" not in compounds.find_one({"inchikey": inchikey}):
            raise InChiKeyNotFoundError(
                "No SMILES feild present for the given InChiKey in the compounds collection of the database"
            )

        if collection.find_one({"inchikey": inchikey}) is not None:
            entry = rdkit.__calc_properties(inchikey, compounds, cipher_mid)
            collection.find_one_and_update({"inchikey": inchikey}, {"$set": {"rdkit": entry}})
            collection.find_one_and_update(
                {"inchikey": inchikey}, {"$currentDate": {"rdkit.modified": True}}
            )
        else:
            data = rdkit.__calc_properties(inchikey, compounds, cipher_mid)
            entry = {"inchikey": inchikey, "rdkit": data}
            db_entry = collection.insert_one(entry)
            collection.find_one_and_update(
                {"inchikey": inchikey}, {"$currentDate": {"rdkit.modified": True}}
            )
            print(db_entry.inserted_id)

    @staticmethod
    def remove():
        """ """
        pass
