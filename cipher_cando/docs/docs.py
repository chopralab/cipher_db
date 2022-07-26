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


def check_cando_id_in_cando(var):
    if Cando.objects.with_id(var) is None:
        raise me.ValidationError("CANDO ID not registered in CANDO collection")


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


def check_sig_id_in_biosig(var):
    if Biosignatures.objects.with_id(var) is None:
        raise me.ValidationError("Biosig ID not registered in Biosignatures collection")


class Biosignatures(me.Document):
    # Biosig ID
    cipher_sig_id = me.StringField(required=True, primary_key=True)
    cipher_mid = me.StringField(required=True, validation=check_mid_in_models)
    description = me.StringField(default="")
    receptors = me.ListField(me.StringField(), required=True, validation=validate_receptor_list)
    scores = me.ListField(me.DecimalField(), required=True)
    modified = me.DateTimeField(default=datetime.datetime.utcnow)


class knn_tuple(me.EmbeddedDocument):
    inchikey = me.StringField(required=True, validation=check_inchikey_in_compounds)
    smiles = me.StringField(required=True, validation=validate_smiles)
    # We can set min and max values here if needed
    cosine_dist = me.DecimalField(required=True)
    rank = me.IntField(required=True)


def check_knn_id_in_knn(var):
    if Knn.objects.with_id(var) is None:
        raise me.ValidationError("KNN ID not registered in biosignatures collection")


class Knn(me.Document):
    cipher_knn_id = me.StringField(required=True, primary_key=True)
    description = me.StringField(default="")
    # Should be set to inchikey if for each compound
    #inchikey = me.StringField(required=True, primary_key=True, validation=check_inchikey_in_compounds)
    #smiles = me.StringField(required=True, validation=validate_smiles)
    cipher_sig_id = me.StringField(required=True, validation=check_sig_id_in_biosig)
    #biosig_scores = me.ListField(me.DecimalField(), required=True)
    # We can probably define the max length of the list if needed
    neighbors = me.ListField(me.EmbeddedDocumentField(knn_tuple), required=True)
    modified = me.DateTimeField(default=datetime.datetime.utcnow)


def gen_unique_cando_id():
    su = shortuuid.ShortUUID(alphabet="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    assigned = False
    while not assigned:
        rand_id = su.random(length=6)
        try: 
            Cando.objects.with_id(rand_id).count()
        except:
            return rand_id


def gen_unique_biosig_id():
    su = shortuuid.ShortUUID(alphabet="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    assigned = False
    while not assigned:
        rand_id = su.random(length=6)
        try:
            Biosignatures.objects.with_id(rand_id).count()
        except:
            return rand_id


def gen_unique_knn_id():
    su = shortuuid.ShortUUID(alphabet="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    assigned = False
    while not assigned:
        rand_id = su.random(length=6)
        try:
            Knn.objects.with_id(rand_id).count()
        except:
            return rand_id
