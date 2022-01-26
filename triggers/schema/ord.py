import pymongo
import json
import requests

class InvalidReactionIDError(Exception):
    pass

class ReactionNotFoundError(Exception):
    pass

class CompoundNotFoundError(Exception):
    pass

class InvalidArgumentError(Exception):
    pass

# TODO: May need additional methods as I am not sure how the ORD fully functions with searching for reactions

def insert(reaction_id):
    '''
    Inserts a specific ORD reaction to the ord collection of the database

    Parameters
    ----------
    reaction_id: string, required
        The ORD reaciton ID of the reaction to be inserted

    Raises
    ------
    InvalidReactionIDError
        TODO: If there is an easy way to check for validity we can throw this error
    '''
    pass

def remove(reaction_id):
    '''
    Removes a specific ORD reaction from the ord collection of the database

    Parameters
    ----------
    reaction_id: string, required
        The ORD reaction ID of the reaction to be removed

    Raises
    ------
    ReactionNotFoundError
        If the ORD reaciton ID is not found in the ord collection of the database
    '''
    pass

def find_all_reactions_for_compound(inchikey,product=True,reactant=False):
    '''
    Returns a list of reaction id of which the compound was involved as either a reactant and/or a product

    Parameters
    ----------
    inchikey: string, required
        The inchikey of the compound in question
    product: boolean, default - True
        True if the search looks for the compound as a product of a reaction
    reactant: boolean, default - False
        True if the search looks for the compound as a reactant of a reaction

    Raises
    ------
    CompoundNotFoundError
        If the compound InChi Key is not found in any of the reactions

    Returns
    -------
    reaction_list: list(string)
        A list of ORD reaction id's of reactions in which the specified compound is involved
    '''
    pass