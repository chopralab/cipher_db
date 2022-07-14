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
    task_id = me.StringField()


class Difficulty(me.Document):
    inchikey = me.StringField(primary_key=True)
    smiles = me.StringField()
    sc_score = me.FloatField()
    sa_score = me.FloatField()
