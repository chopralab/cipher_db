import mongoengine as me
import datetime
from rdkit import Chem


def validate_smiles(val):
    try:
        Chem.MolFromSmiles(val)
    except:
        raise me.ValidationError("Invalid SMILES")


class Compounds(me.Document):
    name = me.StringField(default="")
    id = me.StringField(required=True, primary_key=True)
    inchikey = me.StringField(required=True)
    smiles = me.StringField(required=True, validation=validate_smiles)
    inchi = me.StringField(required=True)
    cid = me.IntField(default=-1)
    iupac = me.StringField(default="")
    modified = me.DateTimeField(default=datetime.datetime.utcnow)


def check_inchikey_in_compounds(val):
    if Compounds.objects.with_id(val) is None:
        raise me.ValidationError("InChI Key not registered in compounds collection")


class Substances(me.Document):
    # Cipher SID
    id = me.StringField(required=True, primary_key=True)
    smiles = me.StringField(required=True, validation=validate_smiles)
    inchikey = me.StringField(required=True, validation=check_inchikey_in_compounds)
    source = me.StringField(required=True)
    name = me.StringField(default="")
    modified = me.DateTimeField(default=datetime.datetime.utcnow)


class Models(me.Document):
    # Cipher MID
    id = me.StringField(required=True, primary_key=True)
    source = me.StringField(required=True)
    parameters = me.StringField(default="\{\}")
    modified = me.DateTimeField(default=datetime.datetime.utcnow)


def check_mid_in_models(val):
    if Models.objects.with_id(val) is None:
        raise me.ValidationError("MID not registered in Models collection")


class External_Databases(me.Document):
    # Cipher EXDBID
    id = me.StringField(required=True, primary_key=True)
    source = me.StringField(required=True)
    url = me.StringField(default="")
    modified = me.DateTimeField(default=datetime.datetime.utcnow)


class Biomolecules(me.Document):
    # Cipher BMID
    id = me.StringField(required=True, primary_key=True)
    name = me.StringField(default="")
    pdb_id = me.StringField(default="")
    chain_id = me.StringField(default="")
    uniprot_id = me.StringField(default="")
    sequence = me.StringField(required=True)
    # TODO: Check to ensure it is in the database
    cipher_mid = me.StringField(required=True, validation=check_mid_in_models)
    modified = me.DateTimeField(default=datetime.datetime.utcnow)


def check_bmid_in_biomolecules(val):
    if Biomolecules.objects.with_id(val) is None:
        raise me.ValidationError("BMID not registered in biomolecules collection")


class Binding_Sites(me.Document):
    # Chiper BSID
    id = me.StringField(required=True, primary_key=True)
    cipher_bmid = me.StringField(required=True, validation=check_bmid_in_biomolecules)
    template_pdb_id = me.StringField(default="")
    template_chain_id = me.StringField(default="")
    residue = me.StringField(default="")
    # TODO: Is this coming in from a file we may need to modify this to a file field
    coordinates = me.StringField(default="")
    modified = me.DateTimeField(default=datetime.datetime.utcnow)


def check_bsid_in_binding_sites(val):
    if Binding_Sites.objects.with_id(val) is None:
        raise me.ValidationError("BSID not registered in binding sites collection")


class Ligands(me.Document):
    # Cipher LID
    id = me.StringField(required=True, primary_key=True)
    inchikey = me.StringField(required=True, validation=check_inchikey_in_compounds)
    pdb_ligand_id = me.StringField(default="")
    # TODO: check this to ensure that it is in the database
    cipher_bsid = me.StringField(required=True, validation=check_bsid_in_binding_sites)
    modified = me.DateTimeField(default=datetime.datetime.utcnow)
