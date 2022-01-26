import pymongo
import json
import requests

class InvalidJsonError(Exception):
    pass

class StructNotFoundError(Exception):
    pass

class InChiKeyNotFoundError(Exception):
    pass

class pubchem:
    '''
    Class for interacting with the PubChem struct of the properties collection of the database

    Methods
    -------
    insert(inchikey)
        Mines physical/chemical property information on the specified compound from the PubChem database and inserts it into the properties collection of the database
    remove(inchikey)
        Removes the PubChem struct of the specified compound from the properties collection of the database
    edit(inchikey, json)
        Edits the PubChem struct of the specified compound based on the JSON information provided
    '''
    @staticmethod
    def insert(inchikey):
        '''
        Mines physical/chemical property information on the specified compound from the PubChem database and inserts it into the properties collection of the database

        Parameters
        ----------
        inchikey: string, required
            The InChi Key of the compound to be mined

        Raises
        ------
        InChiKeyNotFoundError
            If the InChiKey does not exist in the properties collection of the database
        StructNotFoundError
            If a PubChem struct is not found for the given InChi Key entry in the properties collection
        '''
        pass

    @staticmethod
    def remove(inchikey):
        '''
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
        '''
        pass

    @staticmethod
    def edit(inchikey, json):
        '''
        Edits the PubChem struct of the specified compound based on the JSON information provided\
        
        Parameters
        ----------
        inchikey: string, required
            The InChi Key of the entry which has its PubChem struct edited
        json: dict, required
            The JSON doccument which contains edits

        Raises
        ------
        InChiKeyNotFoundError
            If the InChiKey does not exist in the properties collection of the database
        StructNotFoundError
            If a PubChem struct is not found for the given InChi Key entry in the properties collection
        InvalidJsonError
            Conditions yet to be determined
        '''
        pass