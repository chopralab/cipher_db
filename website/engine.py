# Import statements
from decimal import Decimal
from email.mime import base
from email.policy import default
from flask import Flask
from flask_pymongo import PyMongo
from flask_mongoengine import MongoEngine
import mongoengine as me
import datetime
from gridfs import GridFS
import codecs
import base64

#------------------------------------------------------
#--------------- Flask Initilization ------------------
#------------------------------------------------------

# Uncomment this for logging
# logging.basicConfig(filename='logs/restful.log', level=logging.DEBUG)

# Initalize the flask application
app = Flask(__name__)

# Set up Flask-Mongo DB connections
login = open('../utils/login.txt', 'r')
username = login.readline().replace('\n','')
password = login.readline().replace('\n','')
login.close()
mongo_login = 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
app.config['MONGO_URI'] = mongo_login
app.config['MONGODB_SETTINGS'] = {
    'host': 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority',
    'connect': False
}

# Initilize PyMongo
client = PyMongo(app)

# Initilize Mongo Engine
db = MongoEngine()
db.init_app(app)

'''
Compounds Schema
'''

class Compounds(me.Document):
    name = me.StringField(default="")
    inchikey = me.StringField(required=True)
    smiles = me.StringField(required=True)
    inchi = me.StringField(required=True)
    cid = me.IntField(default=-1)
    iupac = me.StringField(default="")
    modified = me.DateTimeField(default=datetime.datetime.utcnow)

'''
Properties Schema
'''

class p_pubchem(me.EmbeddedDocument):
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

class rdkit(me.EmbeddedDocument):
    cipher_mid = me.StringField()
    MolWt = me.DecimalField()
    ExactMolWt = me.DecimalField()
    HeavyAtomMolWt = me.DecimalField()
    MaxPartialCharge = me.DecimalField()
    MinPartialCharge = me.DecimalField()
    NumRadicalElectrons = me.IntField()
    MolLogP = me.DecimalField()
    MaxQED = me.DecimalField()
    MeanQED = me.DecimalField()
    NoneQED = me.DecimalField()
    NHOHCount = me.IntField()
    NOCount = me.IntField()
    NumHAcceptors = me.IntField()
    NumHDonors = me.IntField()
    NumHeteroatoms = me.IntField()
    NumRotateableBonds = me.IntField()
    RingCount = me.IntField()
    modified = me.DateTimeField(default=datetime.datetime.utcnow)

class Properties(me.Document):
    inchikey = me.StringField()
    pubchem = me.EmbeddedDocumentField(p_pubchem)
    rdkit = me.EmbeddedDocumentField(rdkit)

'''
Reactivity Schema
'''

class synthscores(me.EmbeddedDocument):
    sascore = me.DecimalField()
    scscore = me.DecimalField()

class Reactivity(me.Document):
    smiles = me.StringField()
    synthscores = me.EmbeddedDocumentField(synthscores)


'''
Binding Schema
'''

'''
Pubchem Binding Schema
'''

class pubchem_chembl(me.EmbeddedDocument):
    mecid = me.StringField()
    cid = me.StringField()
    chemblid = me.StringField()
    drugname = me.StringField()
    moa = me.StringField()
    action = me.StringField()
    targetchemblid = me.StringField()
    targetname = me.StringField()
    protacxns = me.ListField(me.StringField())
    geneids = me.ListField(me.IntField())
    cmpdname = me.StringField()

class pubchem_drugbank(me.EmbeddedDocument):
    protacxn = me.StringField()
    geneid = me.StringField()
    genesymbol = me.StringField()
    cid = me.StringField()
    drugtype = me.StringField()
    drugname = me.StringField()
    druggroup = me.StringField()
    drugaction = me.StringField()
    drugdetail = me.StringField()
    targettype = me.StringField()
    targetid = me.StringField()
    targetname = me.StringField()
    targetcomponent = me.StringField()
    targetcomponentname = me.StringField()
    generalfunc = me.StringField()
    specificfunc = me.StringField()
    pmids = me.ListField(me.IntField())
    cmpdname = me.StringField()
    dois = me.ListField(me.StringField())

class pubchem_iuphar(me.EmbeddedDocument):
    protacxn = me.StringField()
    geneid = me.StringField()
    ligand = me.StringField()
    cid = me.StringField()
    primarytarget = me.StringField()
    type = me.StringField()
    action = me.StringField()
    units = me.StringField()
    affinity = me.StringField()
    pmids = me.ListField(me.IntField())
    targetid = me.StringField()
    targetname = me.StringField()
    targetspecies = me.StringField()
    genesymbol = me.StringField()
    cmpdname = me.StringField()
    dois = me.ListField(me.StringField())

class pubchem_receptor(me.EmbeddedDocument):
    chembl = me.EmbeddedDocumentField(pubchem_chembl)
    drugbank = me.EmbeddedDocumentField(pubchem_drugbank)
    iuphar = me.EmbeddedDocumentField(pubchem_iuphar)

class b_pubchem(me.EmbeddedDocument):
    mu = me.EmbeddedDocumentField(pubchem_receptor)
    delta = me.EmbeddedDocumentField(pubchem_receptor)
    kappa = me.EmbeddedDocumentField(pubchem_receptor)
    sigma = me.EmbeddedDocumentField(pubchem_receptor)
    nociceptin = me.EmbeddedDocumentField(pubchem_receptor)

'''
DESI Binding Schema
'''

class desi_assay(me.EmbeddedDocument):
    cipher_eid = me.StringField()
    units = me.StringField()
    compound_concentration = me.DecimalField()
    response = me.DecimalField()
    receptor_concentration = me.DecimalField()
    competitor = me.StringField()
    competitor_concentration = me.DecimalField()
    control = me.DecimalField()
    ratio = me.DecimalField()

class desi(me.EmbeddedDocument):
    Mu = me.ListField(desi_assay)
    Delta = me.ListField(desi_assay)
    Sigma = me.ListField(desi_assay)
    Kappa = me.ListField(desi_assay)
    Nociceptin = me.ListField(desi_assay)
    modified = me.DateTimeField(default=datetime.datetime.utcnow)

'''
CANDO Binding Schema
'''

class cando_target(me.EmbeddedDocument):
    cipher_bmid = me.StringField()
    cipher_bsid = me.StringField()
    score = me.DecimalField()

class cando(me.EmbeddedDocument):
    cipher_mid = me.StringField()
    P14416 = me.EmbeddedDocumentField(cando_target)
    P35372 = me.EmbeddedDocumentField(cando_target)
    P35462 = me.EmbeddedDocumentField(cando_target)
    P41143 = me.EmbeddedDocumentField(cando_target)
    P41145 = me.EmbeddedDocumentField(cando_target)
    P41146 = me.EmbeddedDocumentField(cando_target)
    P42262 = me.EmbeddedDocumentField(cando_target)
    Q13224 = me.EmbeddedDocumentField(cando_target)

'''
KNN Binding Schema
'''

class KNN(me.EmbeddedDocument):
    pass

class Binding(me.Document):
    inchikey = me.StringField()
    pubchem = me.EmbeddedDocumentField(b_pubchem)
    desi = me.EmbeddedDocumentField(desi)
    CANDO = me.EmbeddedDocumentField(cando)
    modified = me.DateTimeField(default=datetime.datetime.utcnow)

def return_compounds(identifier):
    '''
    Returns chemical compound information for a compound specified by a specific identifier

    Parameters:
    -----------
    identifier: string, required
        identifer pertaining to the compound in question (Name, SMILES, InChI Key, Pubchem CID, InChI, IUPAC Name)

    Returns:
    --------
    A JSON with all identifying information on the compound loacted in the CIPHER Database or None if it is not found
    '''
    if Compounds.objects(name=identifier).count() > 0:
        return Compounds.objects(name=identifier).to_json()
    elif Compounds.objects(smiles=identifier).count() > 0:
        return Compounds.objects(smiles=identifier).to_json()
    elif Compounds.objects(inchikey=identifier).count() > 0:
        return Compounds.objects(inchikey=identifier).to_json()
    elif Compounds.objects(cid=identifier).count() > 0:
        return Compounds.objects(cid=identifier).to_json()
    elif Compounds.objects(inchi=identifier).count() > 0:
        return Compounds.objects(inchi=identifier).to_json()
    elif Compounds.objects(iupac=identifier).count() > 0:
        return Compounds.objects(iupac=identifier).to_json()
    else:
        return None

def return_properties(inchikey):
    '''
    Returns chemical property information on the compound specified by the given InChI Key

    Parameters:
    -----------
    inchikey: string, required
        The InChI Key of the selected compound

    Returns:
    A JSON with the chemical property information for the selected compound store in the CIPHER database
    '''
    if Properties.objects(inchikey=inchikey).count() > 0:
        return Properties.objects(inchikey=inchikey).to_json()
    else:
        return None

def return_biosignature(inchikey):
    '''
    Return the CANDO biosignature for the selected compound

    Parameters:
    -----------
    inchikey: string, required
        The InChi Key of the selected compound

    Returns:
    The CANDO biosignature of the selected compound in JSON format or None if not found
    '''
    if Binding.objects(inchikey=inchikey).count() > 0:
        return Binding.objects(inchikey=inchikey)[0].CANDO.to_json()
    else:
        return None
        

def return_binding_assays(inchikey):
    '''
    Returns binding assay information on the selected compound

    Parameters:
    -----------
    inchikey: string, required
        The InChI Key of the selected compound

    Returns:
    --------
    pubchem and desi binding assay information in JSON format or None if not found
    '''
    if Binding.objects(inchikey=inchikey).count() > 0:
        pubchem = Binding.objects(inchikey=inchikey)[0].pubchem.to_json()
        desi = Binding.objects(inchikey=inchikey)[0].desi.to_json()
        return pubchem, desi
    else:
        return None, None
    
def return_askcos_pathways(inchikey):
    '''
    Returns all ASKCOS pathways that were found for the given molecule

    Parameters:
    -----------
    inchikey: string, required
        The InChi Key of the selected compound

    Returns:
    -------
    images: list(string)
        A list of base64 decoded image strings, see notes for how to convert the string to an image
    
    Notes:
    ------
    To decode the image use the following:
        decoded_image = open("path1.png", 'wb')
        decoded_image.write(base64.b64decode(images[0]))
        decoded_image.close()
    '''
    fs = GridFS(client.db)
    askcos = client.db.reactivity.find_one({"inchikey": inchikey})['askcos']
    images = []
    for num in askcos['images']:
        image = askcos['images'][num]['imageID']
        gOut = fs.get(image)
        base64_data = codecs.encode(gOut.read(), 'base64')
        image = base64_data.decode('utf-8')
        print(image)
        images.append(image)
    return images

def testing():
    pass

if __name__ == "__main__":
    testing()