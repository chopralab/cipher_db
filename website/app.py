# Import statements
from flask import Flask, request, render_template, redirect, jsonify
from engine import return_compounds, return_properties, return_biosignature, return_askcos_pathways, return_assays
from utils import *

#------------------------------------------------------
#--------------- Flask Initilization ------------------
#------------------------------------------------------

# Uncomment this for logging
# logging.basicConfig(filename='logs/restful.log', level=logging.DEBUG)

# Initalize the flask application
app = Flask(__name__)

# Add this when we deploy to a domain with a SSL certificate
# sslify = SSLify(app)

#---------------------------------------------------
#----------------- Website Backend -----------------
#---------------------------------------------------

# Route to index page
@app.route('/', methods=['GET','POST'])
def index():
    if request.method == 'GET':
        return render_template("index.html"), 200
    elif request.method == 'POST':
        pass
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# Route to the search page
@app.route('/search', methods=['GET','POST'])
def search():
    if request.method == 'GET':
        return render_template('search.html'), 200
    elif request.method == "POST":
        # Replace "naloxone" with the search term from the website frontend search bar
        identifier = request.json['term']
        print(identifier)
        #identifier = "naloxone"
        compounds_id_info = return_compounds(identifier)
        compounds_property_info = []
        compounds_assay_info = []
        compounds_binding_sigs = []
        compounds_retro_pathways = []
        for doc in compounds_id_info:
            print(doc)
            compounds_property_info.append(return_properties(doc["_id"]))
            compounds_assay_info.append(return_assays(doc["_id"]))
            compounds_binding_sigs.append(return_biosignature(doc["_id"]))
            compounds_retro_pathways.append(return_askcos_pathways(doc["_id"]))
            
        return jsonify({"ids": compounds_id_info, "props": compounds_property_info, "biosigs":compounds_binding_sigs, "assays": compounds_assay_info, "synths": compounds_retro_pathways})
        # Compound id info has identifying information about each compound picked up by the search --- list(json formatted dict)
        # Compound property info has propery information from pubchem and rdkit of each compound picked up by the search --- list(json fromatted dict)
        # Compound assay info has assay information from various sources as a JSON formatted dict document --- list(json formatted dict)
        # Compound binding sigs has a list of binding signatures (9 values, none have been assinged in the database yet) for each comound --- list(list(floats))
        # Compound retro pathways has a list of a list of utf-8 encoded image strings corresponding to the ranked retrosynthetic pathways --- list(list(string))
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
        # These are the same as in the search method except one layer of the lists are removed because there is only one compound now
        return render_template("summary.html",name=compounds_id_info[0]["name"].capitalize(),smiles=compounds_id_info[0]["smiles"],inchi=compounds_id_info[0]["inchi"],molformula=render_mol_formula(compounds_property_info["pubchem"]["MolecularFormula"]),molwt=compounds_property_info["pubchem"]["MolecularWeight"],hdc=compounds_property_info["pubchem"]["HBondDonorCount"],hac=compounds_property_info["pubchem"]["HBondAcceptorCount"],logp=compounds_property_info["rdkit"]["MolLogP"]), 200
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

#---------------------------------------------------
#----------------- Login Manager -------------------
#---------------------------------------------------

# Run the Flask Application
if __name__ == "__main__":
    app.run()
