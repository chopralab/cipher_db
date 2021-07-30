# Import statements
from types import MethodType
from flask import Flask, jsonify, request, render_template, redirect
from flask.wrappers import Request
from flask_api import status
from flask_restful import Api
from flask import Flask, request, render_template, redirect
from flask_pymongo import PyMongo
from flask_login import LoginManager, login_required, current_user, login_user, logout_user
from flask_sslify import SSLify
from bson.json_util import dumps
from flask_mongoengine import MongoEngine
import json
import logging
import hashlib

#------------------------------------------------------
#--------------- Flask Initilization ------------------
#------------------------------------------------------

# Uncomment this for logging
# logging.basicConfig(filename='logs/restful.log', level=logging.DEBUG)

# Initalize the flask application
app = Flask(__name__)

# Initilize the Secret Key
with open('utils/secret_key.txt', 'r') as fid:
    app.secret_key = fid.readline().replace('\n','')

# Set up Flask-Mongo DB connections
login = open('utils/login.txt', 'r')
username = login.readline().replace('\n','')
password = login.readline().replace('\n','')
login.close()
mongo_login = 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
app.config['MONGO_URI'] = mongo_login
app.config['MONGODB_SETTINGS'] = {
    'host': 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
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
        return render_template("index.html"), 200
    elif request.method == 'POST':
        if 'gui' in request.form:
            return redirect('/gui/')
        elif 'api' in request.form:
            return redirect('/rest/')
        elif 'login' in request.form:
            return redirect('/login/')
        else:
            return "<pre>" + "Form does not contain correct form info" + "</pre>", 400
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# Route to the login page
@app.route('/login/', methods=['GET','POST'])
def login():
    if request.method == 'GET':
        return render_template("login.html"), 200
    elif request.method == 'POST':
        username = request.form['username']
        password = request.form['password']
        salt_file = open('utils/hash_salt.txt', 'r')
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

# ------------------------------------------------------
# ---------------------- GUI ---------------------------
# ------------------------------------------------------

# Route to the GUI selection pannel
@app.route('/gui/', methods=['GET','POST'])
@login_required
def gui_select():
    if request.method == 'GET':
        return render_template('gui_select.html'), 200
    elif request.method == 'POST':
        if "search" in request.form:
            return redirect('/gui/search/')
        elif "add" in request.form:
            return redirect('/gui/add/')
        elif "update" in request.form:
            return redirect("/gui/update/")
        else:
            return redirect("/")
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# Route to the GUI search engine
@app.route('/gui/search/', methods=['GET','POST'], defaults={'info':None})
@app.route('/gui/search/<info>/', methods=['GET','POST'])
@login_required
def gui_search(info):
    if request.method == 'GET':
        # We can load this as a django webapp later on
        return render_template('gui_search.html'), 200
    elif request.method == 'POST':
        form = request.form
        if "home" in form:
            return redirect("/")
        collection = form['coll']
        input_type = form['input']
        input_value = form['value']
        output_type = form['output']
        params = "?"+"coll="+collection+"&input_type="+input_type+"&input_value="+input_value+"&output_type="+output_type
        return redirect('/gui/search/'+params)
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# Route to the GUI pannel for adding information to the database
@app.route('/gui/add/', methods=['GET','POST'], defaults={'info':None})
@app.route('/gui/add/<info>/', methods=['GET','POST'])
@login_required
def gui_add(info):
    if request.method == 'GET':
        return render_template('gui_add.html')
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# Route to the GUI pannel for updating information in the database
@app.route('/gui/update/', methods=['GET','POST'], defaults={'info':None})
@app.route('/gui/update/<info>/', methods=['GET','POST'])
@login_required
def gui_update(info):
    if request.method == 'GET':
        return render_template('gui_update.html')
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# Route to GUI information pannel
@app.route('/gui/info/')
@login_required
def gui_info():
    if request.method == 'GET':
        return render_template('gui_info.html'), 200
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# -------------------------------------------------
# ---- Query construction for the RESTful API -----
# -------------------------------------------------

# Keywords
# all --> refers to all dataset information for a specific compound (used with dataset)
# summary --> refers to all enteries for the selected dataset (used with data_type)

# Construct a Mongo DB query to access specified information from the database
def construct_query(collection,input_type,input):
    query = {}

    if input_type in ['inchi','inchikey','smiles','name'] and collection == 'general':
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
@login_required
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
        print("Data: " + str(request.form))
        data = request.form
        if "data" not in data:
            return "<pre>" + "Must have \"data\" as a field in your DELETE request" + "</pre>", 400
        try:
            data = json.loads(data['data'])
        except:
            return "<pre>" + "The \"data\" values must be in JSON format" + "</pre>", 400
        doc = client.db.general.delete_one(data)
        print(type(doc))
        return str(data), 200
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# RESTful API URL structure for access to the properties data collection
@app.route('/rest/properties/', methods=['POST','PUT','DELETE'], defaults={'input_type':None,'input':None,'dataset':None,'data_type':None,'output_type':None})
@app.route('/rest/properties/<input_type>/<input>/<dataset>/<data_type>/<output_type>', methods=['GET'])
@login_required
def properties(input_type,input,dataset,data_type,output_type):
    if request.method == 'GET':
        query = construct_query('general',input_type,input)
        projection = construct_projection(collection='properties',dataset=dataset)
        print(query)
        print(projection)
        cursor = client.db.properties.find(query,projection)
        output = format_output(cursor,output_type)
        return output
    elif request.method == 'POST':
        pass
    elif request.method == 'PUT':
        pass
    elif request.method == 'DELETE':
        pass
    else:
        pass

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
        pass
    elif request.method == 'DELETE':
        pass
    else:
        pass

# RESTful API URL structure for access to the reactivity data collection
@app.route('/rest/reactivity/', methods=['POST','PUT','DELETE'], defaults={'input_type':None,'input':None,'dataset':None,'data_type':None,'output_type':None})
@app.route('/rest/reactivity/<input_type>/<input>/<dataset>/<data_type>/<output_type>')
@login_required
def reactivity(input_type,input,dataset,data_type,output_type):
    if request.method == 'GET':
        pass
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
        print("Data: " + str(request.form))
        data = request.form
        if "data" not in data:
            return "<pre>" + "Must have \"data\" as a field in your PUT request" + "</pre>", 400
        try:
            data = json.loads(data['data'])
        except:
            return "<pre>" + "The \"data\" values must be in JSON format" + "</pre>", 400
        doc = client.db.reactivity.insert_one(data)
        print(type(doc))
        return str(data), 200
    elif request.method == 'DELETE':
        print("Data: " + str(request.form))
        data = request.form
        if "data" not in data:
            return "<pre>" + "Must have \"data\" as a field in your DELETE request" + "</pre>", 400
        try:
            data = json.loads(data['data'])
        except:
            return "<pre>" + "The \"data\" values must be in JSON format" + "</pre>", 400
        doc = client.db.reactivity.delete_one(data)
        print(type(doc))
        return str(data), 200
    else:
        return "<pre>" + "Request method not supported" + "</pre>", 400

# RESTful API URL structure for access to the spectra data collection
@app.route('/rest/spectra/', methods=['POST','PUT','DELETE'], defaults={'input_type':None,'input':None,'dataset':None,'data_type':None,'output_type':None})
@app.route('/rest/spectra/<input_type>/<input>/<dataset>/<data_type>/<output_type>')
def spectra(input_type,input,dataset,data_type,output_type):
    if request.method == 'GET':
        pass
    elif request.method == 'POST':
        pass
    elif request.method == 'PUT':
        pass
    elif request.method == 'DELETE':
        pass
    else:
        pass

# RESTful API URL structure for access to the clinical data collection
@app.route('/rest/clinical/', methods=['POST','PUT','DELETE'], defaults={'input_type':None,'input':None,'dataset':None,'data_type':None,'output_type':None})
@app.route('/rest/clinical/<input_type>/<input>/<dataset>/<data_type>/<output_type>')
def clinical(input_type,input,dataset,data_type,output_type):
    if request.method == 'GET':
        pass
    elif request.method == 'POST':
        pass
    elif request.method == 'PUT':
        pass
    elif request.method == 'DELETE':
        pass
    else:
        pass

# RESTful API URL structure for access to the animal data collection
@app.route('/rest/animal/', methods=['POST','PUT','DELETE'], defaults={'input_type':None,'input':None,'dataset':None,'data_type':None,'output_type':None})
@app.route('/rest/animal/<input_type>/<input>/<dataset>/<data_type>/<output_type>')
def animal(input_type,input,dataset,data_type,output_type):
    if request.method == 'GET':
        pass
    elif request.method == 'POST':
        pass
    elif request.method == 'PUT':
        pass
    elif request.method == 'DELETE':
        pass
    else:
        pass

# RESTful API URL structure for access to the machine learning (ML) data collection
@app.route('/rest/ml/', methods=['POST','PUT','DELETE'], defaults={'input_type':None,'input':None,'dataset':None,'data_type':None,'output_type':None})
@app.route('/rest/ml/<input_type>/<input>/<dataset>/<data_type>/<output_type>')
def ml(input_type,input,dataset,data_type,output_type):
    if request.method == 'GET':
        pass
    elif request.method == 'POST':
        pass
    elif request.method == 'PUT':
        pass
    elif request.method == 'DELETE':
        pass
    else:
        pass

# Run the Flask Application
if __name__ == "__main__":
    app.run()