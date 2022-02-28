import mongoengine as me


class ReactionNode(me.EmbeddedDocument):
    smiles = me.StringField()
    # this is List[ChemicalNode] but circular dependency prevents explicitly naming it
    children = me.ListField(me.GenericEmbeddedDocumentField())
    tforms = me.ListField(me.StringField())
    tsources = me.ListField(me.StringField())
    template_score = me.FloatField()
    plausiblity = me.FloatField()
    rank = me.IntField()
    num_examples = me.IntField()
    necessary_reagent = me.StringField()
    precursor_smiles = me.StringField()
    rms_molwt = me.FloatField()
    num_rings = me.IntField()
    scscore = me.FloatField()
    rxn_id = me.StringField()
    forward_score = me.FloatField()
    class_num = me.IntField()
    class_name = me.StringField()


class ChemicalNode(me.EmbeddedDocument):
    smiles = me.StringField()
    rxn = me.EmbeddedDocumentField(ReactionNode)
    ppg = me.FloatField()
    as_reactant = me.IntField()
    as_product = me.IntField()
    terminal = me.BooleanField()
    chemical_id = me.StringField()


class SyntheticTree(me.EmbeddedDocument):
    depth = me.IntField()
    precursor_cost = me.FloatField()
    score = me.FloatField()
    cluster_id = me.IntField()
    root = me.EmbeddedDocumentField(ChemicalNode)
    image = me.ImageField()


class Retrosynthesis(me.Document):
    inchikey = me.StringField(primary_key=True)
    smiles = me.StringField()
    trees = me.EmbeddedDocumentListField(SyntheticTree)


class Difficulty(me.Document):
    inchikey = me.StringField(primary_key=True)
    smiles = me.StringField()
    sc_score = me.FloatField()
    sa_score = me.FloatField()


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
    """Class for running synthetic accessability (SA) score module and inserting information into
    the reactivity collection of the database

    Methods
    -------
    insert(inchikey)
        Runs the SA score module for the compound with the given InChi Key and inserts the results
        into the reactivity collection of the database
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
    """Class for running the SC score module and inserting the resulting information into the
    reactivity collection of the database

    Methods
    -------
    insert(inchikey)
        Runs the SC score module for the compound with the given InChi Key and inserts the results
        into the reactivity collection of the database
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
