import mongoengine as me
import datetime

from cipher_identifiers.docs.docs import (
    validate_smiles,
    check_inchikey_in_compounds,
    check_mid_in_models,
    check_bmid_in_biomolecules,
    check_bsid_in_binding_sites,
)

from cipher_identifiers.docs.docs import Compounds


class Pubchem(me.EmbeddedDocument):
    MolecularFormula = me.StringField()
    MolecularWeight = me.StringField()
    XLogP = me.DecimalField()
    ExactMass = me.StringField()
    MonoisotopicMass = me.StringField()
    TPSA = me.DecimalField()
    Complexity = me.IntField()
    Charge = me.IntField()
    HBondDonorCount = me.IntField()
    HBondAcceptorCount = me.IntField()
    RotatableBondCount = me.IntField()
    HeavyAtomCount = me.IntField()
    IsotopeAtomCount = me.IntField()
    AtomStereoCount = me.IntField()
    DefinedAtomStereoCount = me.IntField()
    UndefinedAtomStereoCount = me.IntField()
    BondStereoCount = me.IntField()
    DefinedBondStereoCount = me.IntField()
    UndefinedBondStereoCount = me.IntField()
    CovalentUnitCount = me.IntField()
    Volume3D = me.DecimalField()
    XStericQuadrupole3D = me.DecimalField()
    YStericQuadrupole3D = me.DecimalField()
    ZStericQuadrupole3D = me.DecimalField()
    FeatureCount3D = me.IntField()
    FeatureAcceptorCount3D = me.IntField()
    FeatureDonorCount3D = me.IntField()
    FeatureAnionCount3D = me.IntField()
    FeatureCationCount3D = me.IntField()
    FeatureRingCount3D = me.IntField()
    FeatureHydrophobeCount3D = me.IntField()
    ConformerModelRMSD3D = me.DecimalField()
    EffectiveRotorCount3D = me.IntField()
    ConformerCount3D = me.IntField()
    Fingerprint2D = me.StringField()
    modified = me.DateTimeField(default=datetime.datetime.utcnow)


class RDKit(me.EmbeddedDocument):
    cipher_mid = me.StringField(validation=check_mid_in_models)
    MolWt = me.DecimalField()
    ExactMolWt = me.DecimalField()
    HeavyAtomMolWt = me.DecimalField()
    MaxPartialCharge = me.DecimalField()
    MinPartialCharge = me.DecimalField()
    NumRadicalElectrons = me.IntField()
    NumValenceElectrons = me.IntField()
    MolLogP = me.DecimalField()
    MaxQED = me.DecimalField()
    MeanQED = me.DecimalField()
    NoneQED = me.DecimalField()
    NHOHCount = me.IntField()
    NOCount = me.IntField()
    NumHAcceptors = me.IntField()
    NumHDonors = me.IntField()
    NumHeteroatoms = me.IntField()
    NumRotatableBonds = me.IntField()
    RingCount = me.IntField()
    modified = me.DateTimeField(default=datetime.datetime.utcnow)


class Properties(me.Document):
    inchikey = me.StringField(required=True, primary_key=True, validation=check_inchikey_in_compounds)
    smiles = me.StringField(required=True, validation=validate_smiles)
    pubchem = me.EmbeddedDocumentField(Pubchem)
    rdkit = me.EmbeddedDocumentField(RDKit)
