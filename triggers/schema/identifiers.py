import pymongo
import json
import requests
import rdkit

class InvalidJsonError(Exception):
    pass

class CompoundNotFoundError(Exception):
    pass

class SubstanceNotFoundError(Exception):
    pass

class ExperimentNotFoundError(Exception):
    pass

class ModelNotFoundError(Exception):
    pass

class BiomoleculeNotFoundError(Exception):
    pass

class ExternalDatabaseNotFoundError(Exception):
    pass

class BindingSiteNotFoundError(Exception):
    pass

class Compounds:
    '''
    Class for interacting with the compounds collection of the CIPHER database

    Methods
    ------
    insert(json)
        Inserts the compound specified by the JSON doccument into the compounds collection of the database
    remove(inchikey)
        Removes the compound with the following InChi Key from the compounds collection of the database
    '''

    def __request_inchikey():
        '''
        Helper method for getting the compounds InChi Key
        '''
        pass

    def __request_inchi():
        '''
        Helper method for getting the compounds InChi string
        '''
        pass

    def __request_smiles():
        '''
        Helper method for getting the compounds SMILES string
        '''
        pass

    def __request_IUPAC_name():
        '''
        Helper method for getting the compounds IUPAC name
        '''
        pass

    @staticmethod
    def insert(json):
        '''
        Inserts the compound specified by the JSON doccument into the compounds collection of the database

        Parameters
        ----------
        json: dict, required
            A JSON doccument which contains either partial or full identifying information on a chemical compound

        Raises
        ------
        InvalidJsonError
            If there is no "inchikey" or "smiles" keys specified in the JSON doccument, then the JSON is rejected 
        '''
        pass

    @staticmethod
    def remove(inchikey):
        '''
        Removes the compound with the following InChi Key from the compounds collection of the database

        Parameters
        ----------
        inchikey: string, required
            The InChi Key of the compound to be removed

        Raises
        -----
        CompoundNotFoundError
            If there is no compound with that InChi Key in the compounds collection of the database
        '''
        pass

class Substances:
    '''
    Class for interacting with the substances collection of the CIPHER database

    Methods
    -------
    insert(json)
        Inserts the substance specified by the JSON doccument into the substances collection of the database
    remove(cipher_sid)
        Removes the substance with the following CIPHER SID from the substances collection of the database
    '''

    @staticmethod
    def insert(json):
        '''
        Inserts the substance specified by the JSON doccument into the substances collection of the database

        Parameters
        ----------
        json: dict, required
            A JSON doccument which contains either partial or full identifying information and a source on a substance

        Raises
        ------
        InvalidJsonError
            If there is no "inchikey"/"smiles" or "source" keys specified in the JSON doccument, then the JSON is rejected 
        '''
        pass

    @staticmethod
    def remove(cipher_sid):
        '''
        Removes the substance with the following CIPHER SID from the compounds collection of the database

        Parameters
        ----------
        cipher_sid: string, required
            The CIPHER SID of the substance to be removed

        Raises
        ------
        SubstanceNotFoundError
            If there is no substance with that CIPHER SID in the substance collection of the database
        '''
        pass

class Experiments:
    '''
    Class for interacting with the experiments collection of the CIPHER database

    Methods
    -------
    insert(json)
        Inserts the experiment specified by the JSON doccument into the substances collection of the database
    remove(cipher_eid)
        Removes the experiment with the following CIPHER EID from the substances collection of the database
    '''
    @staticmethod
    def insert(json):
        '''
        Inserts the experiment specified by the JSON doccument into the experiments collection of the database

        Parameters
        ----------
        json: dict, required
            A JSON doccument which contains information on an experiment

        Raises
        ------
        InvalidJsonError
            If there is no "source" keys specified in the JSON doccument, then the JSON is rejected 
        '''
        pass

    @staticmethod
    def remove(cipher_eid):
        '''
        Removes the experiment with the following CIPHER EID from the experiments collection of the database

        Parameters
        ----------
        cipher_eid: string, required
            The CIPHER EID of the experiment to be removed

        Raises
        ------
        ExperimentNotFoundError
            If there is no experiment with that CIPHER EID in the experiments collection of the database
        '''
        pass

class Models:
    '''
    Class for interacting with the models collection of the CIPHER database

    Methods
    -------
    insert(json)
        Inserts the model specified by the JSON doccument into the models collection of the database
    remove(cipher_mid)
        Removes the model with the following CIPHER MID from the models collection of the database
    '''
    @staticmethod
    def insert(json):
        '''
        Inserts the model specified by the JSON doccument into the models collection of the database

        Parameters
        ----------
        json: dict, required
            A JSON doccument which contains information on a model

        Raises
        ------
        InvalidJsonError
            If there is no "source" keys specified in the JSON doccument, then the JSON is rejected 
        '''
        pass

    @staticmethod
    def remove(cipher_exdbid):
        '''
        Removes the models with the following CIPHER MID from the model collection of the database

        Parameters
        ----------
        cipher_mid: string, required
            The CIPHER MID of the model to be removed

        Raises
        ------
        ExperimentNotFoundError
            If there is no model with that CIPHER MID in the models collection of the database
        '''
        pass

class Biomolecules:
    '''
    Class for interacting with the biomolecule collection of the CIPHER database

    Methods
    -------
    insert(json)
        Inserts the biomolecule specified by the JSON doccument into the biomolecule collection of the database
    remove(cipher_bmid)
        Removes the biomolecule with the following CIPHER BMID from the biomolecule collection of the database
    '''
    @staticmethod
    def insert(json):
        '''
        Inserts the biomolecule specified by the JSON doccument into the biomolecules collection of the database

        Parameters
        ----------
        json: dict, required
            A JSON doccument which contains information on the biomolecule

        Raises
        ------
        InvalidJsonError
            Conditions yet to be determined
        '''
        pass

    @staticmethod
    def remove(cipher_bmid):
        '''
        Removes the biomolecule with the following CIPHER BMID from the biomolecule collection of the database

        Parameters
        ----------
        cipher_bmid: string, required
            The CIPHER BMID of the biomolecule to be removed

        Raises
        ------
        ExperimentNotFoundError
            If there is no biomolecule with that CIPHER BMID in the biomolecules collection of the database
        '''
        pass

class External_Databases:
    '''
    Class for interacting with the external database collection of the CIPHER database

    Methods
    -------
    insert(json)
        Inserts the external database specified by the JSON doccument into the external databases collection of the database
    remove(cipher_exdbid)
        Removes the external database with the following CIPHER EXDBID from the external databases collection of the database
    '''
    @staticmethod
    def insert(json):
        '''
        Inserts the external database specified by the JSON doccument into the external databases collection of the database

        Parameters
        ----------
        json: dict, required
            A JSON doccument which contains information on the external database

        Raises
        ------
        InvalidJsonError
            If there is no "source" keys specified in the JSON doccument, then the JSON is rejected 
        '''
        pass

    @staticmethod
    def remove(cipher_exdbid):
        '''
        Removes the external database with the following CIPHER EXDBID from the external database collection of the database

        Parameters
        ----------
        cipher_bmid: string, required
            The CIPHER EXDBID of the external database to be removed

        Raises
        ------
        ExperimentNotFoundError
            If there is no external database with that CIPHER EXDBID in the external database collection of the database
        '''
        pass

class Binding_Sites:
    '''
    Class for interacting with the binding sites collection of the CIPHER database

    Methods
    -------
    insert(json)
        Inserts the binding site specified by the JSON doccument into the binding site collection of the database
    remove(cipher_exdbid)
        Removes the binding site with the following CIPHER EXDBID from the binding sites collection of the database
    '''
    @staticmethod
    def insert(json):
        '''
        Inserts the binding site specified by the JSON doccument into the binding sites collection of the database

        Parameters
        ----------
        json: dict, required
            A JSON doccument which contains information on the binding site

        Raises
        ------
        InvalidJsonError
            Conditions yet to be determined
        '''
        pass

    @staticmethod
    def remove(cipher_bsid):
        '''
        Removes the binding site with the following CIPHER BSID from the binding site collection of the database

        Parameters
        ----------
        cipher_bmid: string, required
            The CIPHER BSID of the binding site to be removed

        Raises
        ------
        ExperimentNotFoundError
            If there is no binding site with that CIPHER BSID in the binding site collection of the database
        '''
        pass