# Import statements
from crypt import methods
from re import search
from types import MethodType
from flask import Flask, jsonify, request, render_template, redirect
from flask.wrappers import Request
from numpy import identity
from flask_api import status
from flask_restful import Api
from flask import Flask, request, render_template, redirect
from flask_pymongo import PyMongo
from flask_login import LoginManager, login_required, current_user, login_user, logout_user
from flask_sslify import SSLify
from bson.json_util import dumps
from flask_mongoengine import MongoEngine
from gridfs import GridFS
from engine import return_compounds, return_properties, return_biosignature, return_askcos_pathways, return_binding_assays
import json
import logging
import hashlib
import codecs

#------------------------------------------------------
#--------------- Flask Initilization ------------------
#------------------------------------------------------

# Uncomment this for logging
# logging.basicConfig(filename='logs/restful.log', level=logging.DEBUG)

# Initalize the flask application
app = Flask(__name__)

# Initilize the Secret Key
key_file = open('utils/secret_key.txt', 'r')
app.secret_key = key_file.readline().replace('\n','')
key_file.close()

# Set up Flask-Mongo DB connections
login = open('utils/login.txt', 'r')
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

# Add this when we deploy to a domain with a SSL certificate
# sslify = SSLify(app)

# Initilize the Login Manager
login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = 'login'

# User class for Flask sessions
class User(db.Document):
    username = db.StringField()
    password = db.StringField()
    # Display username
    def to_json(self):
        return {"username": username}
    # User is authenticated
    def is_authenticated(self):
        return True
    # User is active
    def is_active(self):
        return True
    # Get the users ID
    def get_id(self):
        return str(self.id)

# Determine the current logged in user for a session
@login_manager.user_loader
def load_user(user_id):
    return User.objects(id=user_id).first()

#---------------------------------------------------
#----------------- Website Backend -----------------
#---------------------------------------------------

# Route to index page
@app.route('/', methods=['GET','POST'])
def index():
    if request.method == 'GET':
        return render_template("home.html"), 200
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
        identifier = "naloxone"
        compounds_id_info = return_compounds(identifier)
        compounds_property_info = []
        compounds_assay_info = []
        compounds_binding_sigs = []
        compounds_retro_pathways = []
        for doc in compounds_id_info:
            compounds_property_info.append(return_properties(doc["inchikey"]))
            compounds_assay_info.append(return_binding_assays(doc["inchikey"]))
            compounds_binding_sigs.append(return_biosignature(doc["inchikey"]))
            compounds_retro_pathways.append(return_askcos_pathways(doc["inchikey"]))

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
        compounds_assay_info = return_binding_assays(inchikey)
        compounds_binding_sigs = return_biosignature(inchikey)
        compounds_retro_pathways = return_askcos_pathways(inchikey)
        # These are the same as in the search method except one layer of the lists are removed because there is only one compound now
        return render_template("summary.html"), 200
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

# Route to the login page
@app.route('/login/', methods=['GET','POST'])
def login():
    if request.method == 'GET':
        return render_template("login.html"), 200
    elif request.method == 'POST':
        username = request.form['username']
        password = request.form['password']
        salt_file = open('../utils/hash_salt.txt', 'r')
        salt = salt_file.readline().replace('\n','')
        key = hashlib.pbkdf2_hmac('sha256', password.encode('utf-8'), salt.encode('utf-8'), 100000, dklen=128)
        key = str(key)[2:-1]
        user = User.objects(username=username,password=key).first()
        if user:
            login_user(user)
            print("User Login: " + username)
            return redirect("/")
        else:
            print("Login Failed with Username: " + username)
            return redirect("/login/")
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

@app.route('/logout/', methods=['GET','POST'])
@login_required
def logout():
    logout_user()
    return redirect('/')

# Route to the RESTful API instruction page
@app.route('/rest/')
def rest():
    return render_template('rest.html'), 200

# -------------------------------------------------
# ---- Query construction for the RESTful API -----
# -------------------------------------------------

# Keywords
# all --> refers to all dataset information for a specific compound (used with dataset)
# summary --> refers to all enteries for the selected dataset (used with data_type)

# Construct a Mongo DB query to access specified information from the database
def construct_query(collection, input_type, input):
    query = {}

    if input_type in ['inchi','inchikey','smiles','name'] and collection in ['general', 'reactivity']:
        query = {input_type: input}
    elif input_type == 'inchikey':
        query = {input_type: input}
    elif input_type in ['inchi','smiles','name'] and collection != 'general':
        cursor = client.db.general.find({input_type: input},{'_id': False, 'inchikey': True})
        inchikey = cursor[0]['inchikey']
        query = {'inchikey': inchikey}
    else:
        print("Datatype not supported")

    return query

# Construct a Mongo DB projection to display specified information from the database
def construct_projection(collection, dataset='all', data_type='summary'):
    projection = {'_id': False}

    if ',' in dataset:
        dataset = dataset.split(',')

    if ',' in data_type:
        data_type = data_type.split(',')

    if collection == 'general':
        if data_type == 'summary':
            return projection
        else:
            if isinstance(data_type, list):
                for elem in data_type:
                    projection[elem] = True
            else:
                projection[data_type] = True
    elif dataset == 'all':
        return projection
    else:
        if isinstance(dataset, list):
            for elem in dataset:
                projection[elem] = 3
        else:
            projection[dataset] = 3
    
    return projection

# Format the database output for rendering to webpage frontend
def format_output(cursor, output_type):
    if output_type == 'json':
        json = []
        for elem in cursor:
            json.append(elem)
        return "<pre>" + dumps(json, indent=2) + "</pre>"
    elif output_type == 'txt':
        data = ''
        for elem in cursor:
            for item in elem.values():
                data += item
        return "<pre>" + data + "</pre>"
    elif output_type == 'csv':
        data = ''
        for elem in cursor:
            for item in elem.values():
                data += item.replace('\n','') + ','
            data = data[:-1]
            data += '\n'
        return '<pre>' + data + '</pre>'
    else:
        print("Output type not supported ...")
        return "<pre>" + "Output type not supported ..." + "</pre>", 400


# -------------------------------------------
# -------------- RESTful API ----------------
# -------------------------------------------
 
# RESTful API URL structure for access to the general data collection
@app.route('/rest/general/', methods=['POST','PUT','DELETE'], defaults={'input_type':None,'input':None,'data_type':None,'output_type':None})
@app.route('/rest/general/<input_type>/<input>/<data_type>/<output_type>', methods=['GET'])
def general(input_type,input,data_type,output_type):
    if request.method == 'GET':
        query = construct_query('general',input_type,input)
        projection = construct_projection(collection='general',data_type=data_type)
        print(query)
        print(projection)
        cursor = client.db.general.find(query,projection)
        output = format_output(cursor,output_type)
        return output, 200
    elif request.method == 'POST':
        print("Data: " + str(request.form))
        data = request.form
        if "query" not in data or "projection" not in data:
            return "<pre>" + "Must have \"query\" and \"projection\" as fields in your POST request" + "</pre>", 400
        try:
            print(type(data['query']))
            print(data['query'])
            print(data['projection'])
            query = json.loads(data['query'])
            projection = json.loads(data['projection'])
        except:
            return "<pre>" + "The \"query\" and \"projection\" values must be in JSON format" + "</pre>", 400
        for key,value in projection.items():
            if value == 'True':
                projection[key] = True
            if value == 'False':
                projection[key] = False
        cursor = client.db.general.find(query,projection)
        output = format_output(cursor,'json')
        output = output.replace("<pre>","")
        output = output.replace("</pre>","")
        return output, 200
    elif request.method == 'PUT':
        print("Data: " + str(request.form))
        data = request.form
        if "data" not in data:
            return "<pre>" + "Must have \"data\" as a field in your PUT request" + "</pre>", 400
        try:
            data = json.loads(data['data'])
        except:
            return "<pre>" + "The \"data\" values must be in JSON format" + "</pre>", 400
        doc = client.db.general.insert_one(data)
        print(type(doc))
        return str(data), 200
    elif request.method == 'DELETE':
       return "<pre>" + "DELETE requests are not supported via the RESTful API <br> If you beleive inaccurate information has been submitted or you would like to make a retraction, please contact the us at ???" + "</pre>", 400
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# RESTful API URL structure for access to the properties data collection
@app.route('/rest/properties/', methods=['POST','PUT','DELETE'], defaults={'input_type':None,'input':None,'dataset':None,'data_type':None,'output_type':None})
@app.route('/rest/properties/<input_type>/<input>/<dataset>/<data_type>/<output_type>', methods=['GET'])
@login_required
def properties(input_type,input,dataset,data_type,output_type):
    if request.method == 'GET':
        query = construct_query('general', input_type, input)
        projection = construct_projection(
            collection='properties',dataset=dataset
        )
        print(query)
        print(projection)
        cursor = client.db.properties.find(query, projection)
        output = format_output(cursor,output_type)
        return output
    elif request.method == 'POST':
        pass
    elif request.method == 'PUT':
        data = request.json
        doc = client.db.spectra.insert_one(data)
        print(type(doc))
        return str(data), 200
    elif request.method == 'DELETE':
        return "<pre>" + "DELETE requests are not supported via the RESTful API <br> If you beleive inaccurate information has been submitted or you would like to make a retraction, please contact the us at ???" + "</pre>", 400
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# RESTful API URL structure for access to the binding data collection  
@app.route('/rest/binding/', methods=['POST','PUT','DELETE'], defaults={'input_type':None,'input':None,'dataset':None,'data_type':None,'output_type':None})
@app.route('/rest/binding/<input_type>/<input>/<dataset>/<data_type>/<output_type>')
@login_required
def binding(input_type,input,dataset,data_type,output_type):
    if request.method == 'GET':
        pass
    elif request.method == 'POST':
        pass
    elif request.method == 'PUT':
        data = request.json
        doc = client.db.spectra.insert_one(data)
        print(type(doc))
        return str(data), 200
    elif request.method == 'DELETE':
        return "<pre>" + "DELETE requests are not supported via the RESTful API <br> If you beleive inaccurate information has been submitted or you would like to make a retraction, please contact the us at ???" + "</pre>", 400
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# RESTful API URL structure for access to the reactivity data collection
@app.route('/rest/reactivity/', methods=['POST','PUT','DELETE'], defaults={'input_type':None,'input':None,'data_type':None,'output_type':None})
@app.route('/rest/reactivity/<input_type>/<path:input>/<data_type>/<output_type>')
@login_required
def reactivity(input_type,input,data_type,output_type):
    if request.method == 'GET':
        if output_type == 'image':
            fs = GridFS(client.db)
            images = []
            askcos = client.db.reactivity.find_one({input_type: input})['askcos']
            for num in askcos['images']:
                image = askcos['images'][num]['imageID']
                gOut = fs.get(image)
                base64_data = codecs.encode(gOut.read(), 'base64')
                image = base64_data.decode('utf-8')
                images.append(image)
            return render_template('trees.html', images=images)
        query = construct_query('reactivity',input_type,input)
        projection = construct_projection(collection='reactivity', data_type=data_type)
        print(query)
        print(projection)
        cursor = client.db.reactivity.find(query,projection)
        output = format_output(cursor,output_type)
        return output
    elif request.method == 'POST':
        print("Data: " + str(request.form))
        data = request.form
        if "query" not in data or "projection" not in data:
            return "<pre>" + "Must have \"query\" and \"projection\" as fields in your POST request" + "</pre>", 400
        try:
            print(type(data['query']))
            print(data['query'])
            print(data['projection'])
            query = json.loads(data['query'])
            projection = json.loads(data['projection'])
        except:
            return "<pre>" + "The \"query\" and \"projection\" values must be in JSON format" + "</pre>", 400
        for key,value in projection.items():
            if value == 'True':
                projection[key] = True
            if value == 'False':
                projection[key] = False
        cursor = client.db.reactivity.find(query,projection)
        output = format_output(cursor,'json')
        output = output.replace("<pre>","")
        output = output.replace("</pre>","")
        return output, 200
    elif request.method == 'PUT':
        data = request.json
        doc = client.db.spectra.insert_one(data)
        print(type(doc))
        return str(data), 200
    elif request.method == 'DELETE':
        return "<pre>" + "DELETE requests are not supported via the RESTful API <br> If you beleive inaccurate information has been submitted or you would like to make a retraction, please contact the us at ???" + "</pre>", 400
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# RESTful API URL structure for access to the properties data collection
@app.route('/rest/ord/', methods=['POST','PUT','DELETE'], defaults={'input_type':None,'input':None,'dataset':None,'data_type':None,'output_type':None})
@app.route('/rest/ord/<input_type>/<input>/<dataset>/<data_type>/<output_type>', methods=['GET'])
@login_required
def ord(input_type, input, dataset, data_type, output_type):
    if request.method == 'GET':
        if input_type != 'inchikey':
            raise ValueError(
                'ORD only searchable by InChiKey key. '
                f'Got input_type="{input_type}"'
            )
        cursor = client.db.ord.find({"products": input}, {"reaction_id"})
        if output_type == 'website':
            html_text = "<pre>"
            for elem in cursor:
                external_url = "https://client.open-reaction-database.org/id/" + str(elem["reaction_id"])
                reaction_id = str(elem["reaction_id"])
                html_text += "<a href=" + external_url + ">" + reaction_id + "</a>"
                html_text += "<br>"
            html_text +="</pre>"
            return html_text, 200
        else:
            output = format_output(cursor, output_type)
            return output
    elif request.method == 'POST':
        pass
    elif request.method == 'PUT':
        data = request.json
        doc = client.db.spectra.insert_one(data)
        print(type(doc))
        return str(data), 200
    elif request.method == 'DELETE':
         return "<pre>" + "DELETE requests are not supported via the RESTful API <br> If you beleive inaccurate information has been submitted or you would like to make a retraction, please contact the us at ???" + "</pre>", 400
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# RESTful API URL structure for access to the spectra data collection
@app.route('/rest/spectra/', methods=['POST','PUT','DELETE'], defaults={'input_type':None,'input':None,'dataset':None,'data_type':None,'output_type':None})
@app.route('/rest/spectra/<input_type>/<input>/<dataset>/<data_type>/<output_type>')
@login_required
def spectra(input_type,input,dataset,data_type,output_type):
    if request.method == 'GET':
        pass
    elif request.method == 'POST':
        pass
    elif request.method == 'PUT':
        data = request.json
        doc = client.db.spectra.insert_one(data)
        print(type(doc))
        return str(data), 200
    elif request.method == 'DELETE':
         return "<pre>" + "DELETE requests are not supported via the RESTful API <br> If you beleive inaccurate information has been submitted or you would like to make a retraction, please contact the us at ???" + "</pre>", 400
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# RESTful API URL structure for access to the clinical data collection
@app.route('/rest/clinical/', methods=['POST','PUT','DELETE'], defaults={'input_type':None,'input':None,'dataset':None,'data_type':None,'output_type':None})
@app.route('/rest/clinical/<input_type>/<input>/<dataset>/<data_type>/<output_type>')
@login_required
def clinical(input_type,input,dataset,data_type,output_type):
    if request.method == 'GET':
        return "<pre>" + "This database collection has not been implemented yet" + "</pre>", 200
    elif request.method == 'POST':
        return "<pre>" + "This database collection has not been implemented yet" + "</pre>", 200
    elif request.method == 'PUT':
        return  "<pre>" + "Addition to this collection has not been implemented yet" + "</pre>", 400
    elif request.method == 'DELETE':
         return "<pre>" + "DELETE requests are not supported via the RESTful API <br> If you beleive inaccurate information has been submitted or you would like to make a retraction, please contact the us at ???" + "</pre>", 400
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# RESTful API URL structure for access to the animal data collection
@app.route('/rest/animal/', methods=['POST','PUT','DELETE'], defaults={'input_type':None,'input':None,'dataset':None,'data_type':None,'output_type':None})
@app.route('/rest/animal/<input_type>/<input>/<dataset>/<data_type>/<output_type>')
@login_required
def animal(input_type,input,dataset,data_type,output_type):
    if request.method == 'GET':
        return "<pre>" + "This database collection has not been implemented yet" + "</pre>", 200
    elif request.method == 'POST':
        return "<pre>" + "This database collection has not been implemented yet" + "</pre>", 200
    elif request.method == 'PUT':
        return  "<pre>" + "Addition to this collection has not been implemented yet" + "</pre>", 400
    elif request.method == 'DELETE':
         return "<pre>" + "DELETE requests are not supported via the RESTful API <br> If you beleive inaccurate information has been submitted or you would like to make a retraction, please contact the us at ???" + "</pre>", 400
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# RESTful API URL structure for access to the machine learning (ML) data collection
@app.route('/rest/ml/', methods=['POST','PUT','DELETE'], defaults={'input_type':None,'input':None,'dataset':None,'data_type':None,'output_type':None})
@app.route('/rest/ml/<input_type>/<input>/<dataset>/<data_type>/<output_type>')
@login_required
def ml(input_type,input,dataset,data_type,output_type):
    if request.method == 'GET':
        return "<pre>" + "This database collection has not been implemented yet" + "</pre>", 200
    elif request.method == 'POST':
        return "<pre>" + "This database collection has not been implemented yet" + "</pre>", 200
    elif request.method == 'PUT':
        return  "<pre>" + "Addition to this collection has not been implemented yet" + "</pre>", 400
    elif request.method == 'DELETE':
         return "<pre>" + "DELETE requests are not supported via the RESTful API <br> If you beleive inaccurate information has been submitted or you would like to make a retraction, please contact the us at ???" + "</pre>", 400
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# Run the Flask Application
if __name__ == "__main__":
    app.run()