from mongoengine import *

# from triggers.schema.askcos import TreeBuilder, tree_to_image


class ChemicalNode(EmbeddedDocument):
    """empty class so the code can be interpreted"""


class ReactionNode(EmbeddedDocument):
    smiles = StringField()
    tforms = ListField(StringField)
    tsources = ListField(StringField)
    template_score = FloatField()
    plausiblity = FloatField()
    rank = IntField()
    num_examples = IntField()
    necessary_reagent = StringField()
    precursor_smiles = StringField()
    rms_molwt = FloatField()
    num_rings = IntField()
    scscore = FloatField()
    rxn_id = StringField()
    chidren = ListField(EmbeddedDocumentField(ChemicalNode))
    forward_score = FloatField()
    class_num = IntField()
    class_name = StringField()


class ChemicalNode(EmbeddedDocument):
    smiles = StringField()
    ppg = FloatField()
    as_reactant = IntField()
    as_product = IntField()
    terminal = BooleanField()
    chemical_id = StringField()
    children = ListField(EmbeddedDocumentField(ReactionNode))

    meta = {"allow_inheritance": True}


class TreeRoot(ChemicalNode):
    depth = IntField()
    precursor_cost = FloatField()
    score = FloatField()
    cluster_id = IntField()


class Retrosynthesis(Document):
    inchikey = StringField()
    smiles = StringField()
    trees = ListField(EmbeddedDocumentField(TreeRoot))


class Feasibility(Document):
    inchikey = StringField()
    smiles = StringField()
    sc_score = FloatField()
    sa_score = FloatField()


class askcos:
    """
    Class for interacting with ASKCOS structs in the reactivity collection of the database

    Methods
    -------
    insert()
    remove()
    edit()
    """

    @staticmethod
    def insert():
        """ """
        pass

    @staticmethod
    def remove():
        """ """
        pass

    @staticmethod
    def edit():
        """ """
        pass


class sa_score:
    """
    Class for running synthetic accessability (SA) score module and inserting information into the reactivity collection of the database

    Methods
    -------
    insert(inchikey)
        Runs the SA score module for the compound with the given InChi Key and inserts the results into the reactivity collection of the database
    remove(inchikey)
        Remove the SA score struct for the compound with the given InChi Key
    """

    @staticmethod
    def insert(inchikey):
        """ """
        pass

    @staticmethod
    def remove(inchikey):
        """ """
        pass


class sc_score:
    """
    Class for running the SC score module and inserting the resulting information into the reactivity collection of the database

    Methods
    -------
    insert(inchikey)
         Runs the SC score module for the compound with the given InChi Key and inserts the results into the reactivity collection of the database
    remove(inchikey)
        Remove the SC score struct for the compound with the given InChi Key
    """

    @staticmethod
    def insert(inchikey):
        """ """
        pass

    @staticmethod
    def remove(inchikey):
        """ """
        pass
