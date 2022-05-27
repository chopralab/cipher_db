# Import statements
from flask import Flask, request, render_template, redirect, jsonify, session, abort
from engine import return_compounds, return_properties, return_biosignature, return_askcos_pathways, return_assays, return_compound_image, return_desired_dynamic_biosignature, return_askcos_pathways, return_biosig_knn
from utils import *
from json import dumps, loads

import sys
sys.path.append("../")

from cipher_identifiers.utils.compounds import id_compound_from_smiles
from cipher_properties.utils.properties import insert_properties_from_smiles, Properties

#------------------------------------------------------
#--------------- Flask Initilization ------------------
#------------------------------------------------------

# Uncomment this for logging
# logging.basicConfig(filename='logs/restful.log', level=logging.DEBUG)

# Initalize the flask application
app = Flask(__name__)
app.secret_key = 'cipher'

# Add this when we deploy to a domain with a SSL certificate
# sslify = SSLify(app)

#---------------------------------------------------
#----------------- Website Backend -----------------
#---------------------------------------------------

# Route to index page
@app.route('/', methods=['GET','POST'])
def index():
    if request.method == 'GET':
        if "prev_term" not in session: session["prev_term"] = ""
        return render_template("index.html"), 200
    elif request.method == 'POST':
        pass
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400
        
@app.route('/add', methods=['GET','POST'])
def add():
    if request.method == 'GET':
        return render_template("add.html"), 200
    elif request.method == 'POST':
        #use request.json["smiles"] to access the smiles
        request.json["smiles"]
        return jsonify({})
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

@app.route('/top', methods=['GET','POST'])
def top():
    if request.method == 'POST':
        return jsonify({"desired": return_biosig_knn()[0], "top": return_biosig_knn()[1]})
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# Route to the search page
@app.route('/search', methods=['GET','POST'])
def search():
    if request.method == 'GET':
        if "prev_term" not in session: session["prev_term"] = ""
        return render_template('search.html',term=session["prev_term"]), 200
    elif request.method == "POST":
        # Replace "naloxone" with the search term from the website frontend search bar
        identifier = request.json["term"]
        session["prev_term"] = identifier
        print(identifier)
        #identifier = "naloxone"
        compounds_id_info = return_compounds(identifier)
        compounds_property_info = []
        compounds_assay_info = []
        compounds_binding_sigs = []
        compounds_retro_pathways = []
        for doc in compounds_id_info:
            # print(doc)
            compounds_property_info.append(return_properties(doc["_id"]))
            compounds_assay_info.append(return_assays(doc["_id"]))
            compounds_binding_sigs.append(return_biosignature(doc["_id"]))
            compounds_retro_pathways.append(return_askcos_pathways(doc["_id"]))
        for entry in compounds_property_info:
            if "MolecularFormula" in entry["pubchem"]:
                entry["pubchem"]["MolecularFormula"] = render_mol_formula(entry["pubchem"]["MolecularFormula"])
            entry["svg"] = return_compound_image(entry["_id"]).replace("height='300px'","height='200px'").replace("width='300px'","width='200px'")
            #temp appending of pathway images to compounds_retro_pathways:
            compounds_retro_pathways.append(return_askcos_pathways(entry["_id"]))
        respone = {"ids": compounds_id_info, "props": compounds_property_info, "biosigs":compounds_binding_sigs, "assays": compounds_assay_info, "synths": compounds_retro_pathways, "desired": return_desired_dynamic_biosignature()}
        # print(respone)
        return jsonify(respone)
        # Compound id info has identifying information about each compound picked up by the search --- list(json formatted dict)
        # Compound property info has propery information from pubchem and rdkit of each compound picked up by the search --- list(json fromatted dict)
        # Compound assay info has assay information from various sources as a JSON formatted dict document --- list(json formatted dict)
        # Compound binding sigs has a list of binding signatures (9 values, none have been assinged in the database yet) for each comound --- list(list(floats))
        # Compound retro pathways has a list of a list of utf-8 encoded image strings corresponding to the ranked retrosynthetic pathways --- list(list(string))
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

@app.route('/info', methods=['GET','POST'])
def info():
    if request.method == 'POST':
        print(request.json)
        compounds_id_info = return_compounds(request.json['inchikey'])
        compounds_property_info = return_properties(request.json['inchikey'])
        compounds_assay_info = return_assays(request.json['inchikey'])
        compounds_binding_sigs = return_biosignature(request.json['inchikey'])
        compounds_retro_pathways = return_askcos_pathways(request.json['inchikey'])
        if "MolecularFormula" not in compounds_property_info["pubchem"]:
            compounds_property_info["pubchem"]["MolecularFormula"] = ""
        else:
            compounds_property_info["pubchem"]["MolecularFormula"] = render_mol_formula(compounds_property_info["pubchem"]["MolecularFormula"])
        if "MolecularWeight" not in compounds_property_info["pubchem"]:
            compounds_property_info["pubchem"]["MolecularWeight"] = ""
        if "HBondDonorCount" not in compounds_property_info["pubchem"]:
            compounds_property_info["pubchem"]["HBondDonorCount"] = ""
        if "HBondAcceptorCount" not in compounds_property_info["pubchem"]:
            compounds_property_info["pubchem"]["HBondAcceptorCount"] = ""
        compounds_property_info["svg"] = return_compound_image(compounds_property_info["_id"])
        respone = {"ids": compounds_id_info, "props": compounds_property_info, "biosigs":compounds_binding_sigs, "assays": compounds_assay_info, "synths": compounds_retro_pathways, "desired": return_desired_dynamic_biosignature()}
        return jsonify(respone)
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

@app.route('/summary/<inchikey>')
def summary(inchikey):
    if request.method == 'GET':
        compounds_id_info = return_compounds(inchikey)
        compounds_property_info = return_properties(inchikey)
        compounds_assay_info = return_assays(inchikey)
        compounds_binding_sigs = return_biosignature(inchikey)
        compounds_retro_pathways = return_askcos_pathways(inchikey)
        print(compounds_id_info)
        if "MolecularFormula" not in compounds_property_info["pubchem"]:
            compounds_property_info["pubchem"]["MolecularFormula"] = ""
        if "MolecularWeight" not in compounds_property_info["pubchem"]:
            compounds_property_info["pubchem"]["MolecularWeight"] = ""
        if "HBondDonorCount" not in compounds_property_info["pubchem"]:
            compounds_property_info["pubchem"]["HBondDonorCount"] = ""
        if "HBondAcceptorCount" not in compounds_property_info["pubchem"]:
            compounds_property_info["pubchem"]["HBondAcceptorCount"] = ""
        # These are the same as in the search method except one layer of the lists are removed because there is only one compound now
        return render_template("summary.html",name=compounds_id_info[0]["name"].capitalize(),smiles=compounds_id_info[0]["smiles"],inchi=compounds_id_info[0]["inchi"],molformula=render_mol_formula(compounds_property_info["pubchem"]["MolecularFormula"]),molwt=compounds_property_info["pubchem"]["MolecularWeight"],hdc=compounds_property_info["pubchem"]["HBondDonorCount"],hac=compounds_property_info["pubchem"]["HBondAcceptorCount"],logp=compounds_property_info["rdkit"]["MolLogP"],inchikey=inchikey), 200
    elif request.method == 'POST':
        return "<pre>" + "Request method not supported" + "</pre>", 400
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

@app.route('/contact', methods=['GET', 'POST'])
def contact():
    if request.method == 'GET':
        return "<pre>" + "Request method not supported" + "</pre>", 400
    elif request.method == 'POST':
        return "<pre>" + "Request method not supported" + "</pre>", 400
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

@app.route('/rest/properties/<inchikey>', methods=['GET', 'POST'])
def rest_properties(inchikey):
    comp = Properties.objects().with_id(inchikey)
    if comp is not None:
        return "<pre>" + dumps(loads(comp.to_json()), indent=2) + "</pre>"
    else:
        abort(404)
#---------------------------------------------------
#----------------- Login Manager -------------------
#---------------------------------------------------

# Run the Flask Application
if __name__ == "__main__":
    app.run()
