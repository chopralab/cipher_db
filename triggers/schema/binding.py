import pymongo
import json
import requests

class InvalidJsonError(Exception):
    pass

class InvalidFileError(Exception):
    pass

class InChiKeyNotFoundError(Exception):
    pass

class StructNotFoundError(Exception):
    pass

class pubchem():
    '''
    Class for interacting with PubChem binding information in the binding collection of the database

    Methods
    -------
    def insert_from_file()
    def add(inchikey)
    def remove(inchikey)
    def edit(inchikey, json)
    '''
    @staticmethod
    def insert_from_file(fname):
        '''
        '''
        pass

    @staticmethod
    def insert(inchikey):
        '''
        '''
        pass

    @staticmethod
    def remove(inchikey):
        '''
        '''
        pass

    @staticmethod
    def edit(inchikey, json):
        '''
        '''
        pass

class cando():
    '''
    Class for interacting with CANDO binding data in the binding collection of the database

    Methods
    -------
    insert(inchikey)
        Runs the CANDO module on the specified compound and updates the binding collection of the database with the resulting information
    remove(inchikey)
        Removes the CANDO struct of the specified compound from the binding collection of the database
    edit(inchikey, json)
        Edits the CANDO struct of the specified compound with the information in the specified JSON file 
    '''
    @staticmethod
    def insert(inchikey):
        '''
        '''
        pass

    @staticmethod
    def remove(inchikey):
        '''
        '''
        pass

    @staticmethod
    def edit(inchikey, json):
        '''
        '''
        pass

class desi_exp():
    '''
    Class for interacting with DESI experiments in the binding collection of the database

    Methods
    -------
    insert(json)
        Inserts DESI experiment information specified in the JSON file to the binding collection of the database
    remove(inchikey)

    edit(inchikey, json)
    
    '''
    @staticmethod
    def insert(json):
        '''
        '''
        pass

    @staticmethod
    def remove(inchikey):
        '''
        '''
        pass

    @staticmethod
    def edit(inchikey, json):
        '''
        '''
        pass