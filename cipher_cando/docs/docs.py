import mongoengine as me
import sys
import datetime
import shortuuid

sys.path.append("../../")
from cipher_identifiers.docs.docs import (
    validate_smiles,
    check_inchikey_in_compounds,
    check_mid_in_models,
    check_bmid_in_biomolecules,
    check_bsid_in_binding_sites,
)


def validate_receptor_list(var):
    for receptor in var:
        check_bmid_in_biomolecules(receptor)


class Cando(me.Document):
    #Cando ID
    cipher_cando_id = me.StringField(required=True, primary_key=True)
    # Can assign a default MID for CANDO runs if needed, check mongo engine doccumentation
    inchikey = me.StringField(required=True)
    smiles = me.StringField(required=True)
    cipher_mid = me.StringField(required=True, validation=check_mid_in_models)
    cipher_bmid = me.StringField(required=True, validation=check_bmid_in_biomolecules)
    cipher_bsid = me.StringField(required=True, validation=check_bsid_in_binding_sites)
    # Can contrain a min and max value, check mongo engine doccumentation
    interaction_score = me.DecimalField(required=True)
    modified = me.DateTimeField(default=datetime.datetime.utcnow)


class Biosignatures(me.Document):
    # Biosig ID
    cipher_sig_id = me.StringField(required=True, primary_key=True)
    description = me.StringField(default="")
    receptors = me.ListField(me.StringField(), required=True, validation=validate_receptor_list)


def check_sig_id_in_biosig(var):
    if Biosignatures.compounds.with_id(var).count() == 0:
        raise me.ValidationError("Bio Sig ID not registered in biosignatures collection")


class knn_tuple(me.EmbeddedDocument):
    inchikey = me.StringField(required=True, validation=check_inchikey_in_compounds)
    smiles = me.StringField(required=True, validation=validate_smiles)
    # We can set min and max values here if needed
    cosine_dist = me.DecimalField(required=True)


class KNN(me.Document):
    # Should be set to inchikey if for each compound
    inchikey = me.StringField(required=True, primary_key=True, validation=check_inchikey_in_compounds)
    smiles = me.StringField(required=True, validation=validate_smiles)
    cipher_sig_id = me.StringField(required=True, validation=check_sig_id_in_biosig)
    # We can probably define the max length of the list if needed
    neighbors = me.ListField(me.EmbeddedDocumentField(knn_tuple), required=True)

def gen_unique_cando_id():
    su = shortuuid.ShortUUID(alphabet="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    assigned = False
    while not assigned:
        rand_id = su.random(length=6)
        #if CANDO.objects.with_id(rand_id).count() == 0:
        #    return rand_id
        try: 
            CANDO.objects.with_id(rand_id).count()
        except:
            return rand_id

def gen_unique_biosig_id():
    su = shortuuid.ShortUUID(alphabet="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    assigned = False
    while not assigned:
        rand_id = su.random(length=6)
        #if Biosignatures.objects.with_id(rand_id).count() == 0:
        #    return rand_id
        try:
            Biosignatures.objects.with_id(rand_id).count()
        except:
            return rand_id
