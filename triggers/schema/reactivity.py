from tkinter.tix import Tree
from triggers.schema.askcos import TreeBuilder, tree_to_image 


tb = TreeBuilder()

class askcos():
    '''
    Class for interacting with ASKCOS structs in the reactivity collection of the database

    Methods
    -------
    insert()
    remove()
    edit()
    '''
    @staticmethod
    def insert():
        '''
        '''
        pass

    @staticmethod
    def remove():
        '''
        '''
        pass

    @staticmethod
    def edit():
        '''
        '''
        pass

class sa_score():
    '''
    Class for running synthetic accessability (SA) score module and inserting information into the reactivity collection of the database

    Methods
    -------
    insert(inchikey)
        Runs the SA score module for the compound with the given InChi Key and inserts the results into the reactivity collection of the database
    remove(inchikey)
        Remove the SA score struct for the compound with the given InChi Key
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

class sc_score():
    '''
    Class for running the SC score module and inserting the resulting information into the reactivity collection of the database

    Methods
    -------
    insert(inchikey)
         Runs the SC score module for the compound with the given InChi Key and inserts the results into the reactivity collection of the database
    remove(inchikey)
        Remove the SC score struct for the compound with the given InChi Key
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