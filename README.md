# The CIPHER Database

The CIPHER Database is a semi-autonomous information system designed to facilitate the discovery and synthesis of novel, non-addictive opioids. CIPHER is designed to mine and calculate chemical descriptors for novel compounds added to the database in an automated manner, calculate an interaction signature for compounds using the CANDOCK module, attempt to find a retrosynthetic pathway for compounds using the ASKCOS module, and use this information determine an optimum opioid binding profile for both efficacy and non-additivity. CIPHER is also designed to directly interface with any analytical instrument which is connected to the internet and can run a script in any language which can make a POST request to a URL. See the instrumental setup guide on how to connect an instrument to the database. The final major feature of CIPHER is integration of generative models which can use database information to design novel opioid-like compounds based on the current non-addictive opioid profile in the database.

## [Setup](#setup)

The following guide will show you the basics of how to set up the cipher database on a local computer or a cloud hosted machine. Additional instructions are provided for modifying the functionality of the database to account for different drug discovery project. 

### [Requirements](#requirements)

```yaml
Flask==2.0.1
Flask_API==3.0.post1
Flask_Login==0.5.0
flask_mongoengine==1.0.0
Flask_PyMongo==2.3.0
Flask_RESTful==0.3.9
Flask_SSLify==0.1.5
mongoengine==0.23.1
pymongo==3.11.4
```

**Note that this is just a basic list of requirements for general database interaction. To download all of the required packages and set up external modules correctly, please create a conda environment from the `cipher_db_env.yml` file and then navigate over to the specific module page for additional setup instructions**

### Installation
1. `conda env create -f env.yml`
1. `pip install -e .`
1. (optional) install additional requirements for development or deployment: `pip install -e .[OPTION] --upgrade`

possible values for `OPTION`:
* `dev`: for general development purposes (linting, testing, etc.)
* `reactivity`: for reactivity trigger deployment


### [Database Deployment](#database-deployment)

The CIPHER database is currently implemented using the MongoDB database architecture. While it is possible to transfer the data storage mechanism to a different non-relational or relational database architecture, this guide will cover how to set up, maintain, and expand the data storage mechanism in Mongo DB. Please note that there are two potential options for setting up MongoDB, locally hosted and cloud hosted. The MongoDB version which was during development and testing is `v5.0.9` and the storage was cloud hosted, via Microsoft Azure using the M2 (General) Tier. 

#### [Database Setup - Cloud Hosted](#database-setup---cloud-hosted)

To set up MongoDB with a cloud hosted database, please follow the [instructions on how to set up a cloud hosted database using MongoDB Atlas](https://www.mongodb.com/docs/atlas/getting-started/). You will need to create a database account prior to attempting to access the database. Please note that database permissions are handled on MongoDB Atlas, so if you wish to grant certain users certain permissions you would do so here. Additionally, not all features may work properly if the account used does not have the appropriate privileges. The accounts used during development and testing had `readWriteAnyDatabase` permissions with access to all resources.

From here, you will need the obtain the [MongoDB URI](https://www.mongodb.com/docs/atlas/tutorial/connect-to-your-cluster/#connect-to-your-atlas-cluster) to connect to the cloud hosted database on the atlas cluster. After obtaining the MongoDB URI, you will need to set up an environment variable on your local machine which contains the MongoDB application driver with the appropriate username and password. There are two URI's which you can set, one is for a production database and the other is for a testing database. You can create a bash script with the following lines of code to set both URI's easily:

```bash
export MONGO_URI='(URI for Development Database)'
export TESTING_URI='(URI for Testing Database)'
```

From here, your database should be appropriately set up and ready to use.

#### [Database Setup - Locally Hosted](#database-setup---locally-hosted)

To set up MongoDB with a locally hosted database, please navigate over to the [installation tutorial page on the MongoDB website](https://www.mongodb.com/docs/manual/installation/#mongodb-installation-tutorials) and select the appropriate version and operating system (supported options are Windows, MacOS, and various versions of Linux). From here, follow the tutorial instructions to install MongoDB, install mongosh, set up the database directory, and start up the database. [By default, a database which is locally hosted will only be able to accept connections from clients which are run on the locally hosted server](https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-windows/#localhost-binding-by-default). If you wish to allow for a remote client to access the database you will have to change this in the [following manner](https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-windows/#localhost-binding-by-default). It is also recommended that you [install MongoDB Compass](https://www.mongodb.com/try/download/compass) if the local machine has a GUI as this can be used to easily manage the database.

To connect to this database, you will have to obtain a URI connection string for the locally hosted database. How this string is obtained depends on how the database is set up. Please see [this page on MongoDB URI's](https://www.mongodb.com/docs/manual/reference/connection-string/) and view the URI for the deployment type which is used. Other resources such as [MongoDB Compass](https://www.mongodb.com/try/download/compass) may be able to assist in the process of obtaining the database URI.

After obtaining the MongoDB URI, you will need to set up an environment variable on your local machine which contains the MongoDB application driver with the appropriate username and password. There are two URI's which you can set, one is for a production database and the other is for a testing database. You can create a bash script with the following lines of code to set both URI's easily:

```bash
export MONGO_URI='(URI for Development Database)'
export TESTING_URI='(URI for Testing Database)'
```

From here, your database should be appropriately set up and ready to use.

### [Trigger Deployment - Database](#trigger-deployment---database)

The CIPHER Database is designed to automatically generate chemical records which are as complete as possible upon the addition of novel chemical data. This is accomplished using what is known as a MongoDB database trigger. How this works is that when specific chemical identifiers (SMILES, InChI, InChI Key) are added to the database in any entry the following process occurs:

1. A database trigger script checks to see if a general record of the compound exists (database novelty)
    - **General Entry Exists:** No further action is taken except adding the original record which activated the trigger - **END**
    - **General Entry Does Not Exist:** A general entry is added for the compound with basic identifying information (Name (if given), SMILES, InChI Key, Molecular Formula) - **CONTINUE TO STEP 2**
2. A second trigger attempts to mine chemical property information from PubChem (if available) using the REST API and calculatates chemical descriptors using RDKit. This information is stored in the identifies collection of the database
3. A third trigger then attempts to calculate a biosignature using the [CANDO software package](https://github.com/ram-compbio/CANDO). **Note:** The biosignature by default is a series of 8 receptors identified as involved in opioid efficacy and addictively, this can be tuned to a different subset of the human proteome in the following manner (**Link to CANDO README HERE**). 
4. Finally, a trigger script attempts to find a retrosynthetic pathway for the compound using the [ASKCOS software package](https://github.com/ASKCOS/ASKCOS). The software is run using a set of default parameters which can be changed in the following manner (**Link to ASKCOS README HERE**)

To set up a database trigger for use, follow the trigger setup guide for the specific trigger in its own README (links - compounds, properties, assays, CANDO, ASKCOS) and then run the `trigger.py` script from the command line. The script will run indefinitely and may output logging information based on what database operations are encountered and what subsequent actions are preformed. The database trigger modules are designed to run separately and work in parallel. This means that different triggers can be run on separate machines at the same time provide each machine has the appropriate drivers and packages installed to connect to database and run the trigger scripts.

### [Trigger Deployment - Scheduled](#trigger-deployment---scheduled)

While database triggers do an excellent job of keeping the database up to date whenever novel data is added, there are certain fields which only need to be updated periodically. An example of such a field in this project is the desired dynamic biosignature. Updating these fields involves using a MongoDB schedule trigger, which is a segment of code which performs a specific operation on the database after a certain amount of time has passed (e.g. every 24 hours).

## [Modifications](#modifications)

### [Database Collections](#database-collections)

While the main driver for accessing MongoDB database information in Python is PyMongo, this project mainly uses an extension of PyMongo called MongoEngine. MongoEngine allows for database collection documents to be defined via Python classes. More information on how MongoEngine is used can be found [here](https://docs.mongoengine.org/index.html). For example, below is a segment of the `docs.py` file which has the MongoEngine document class for compound identification information.

```python
class Compounds(me.Document):
    inchikey = me.StringField(required=True, primary_key=True)
    name = me.StringField(default="")
    smiles = me.StringField(required=True, validation=validate_smiles)
    inchi = me.StringField(required=True)
    cid = me.StringField(default="")
    iupac = me.StringField(default="")
    synonyms = me.ListField(me.StringField(), default = [])
    image = me.FileField()
    modified = me.DateTimeField(default=datetime.datetime.utcnow)
```

If you wish to record additional information for a specific collection or you wish to remove recorded information, you can add and remove fields from the document. **Note that if you add a required field or remove a field which has its values present in the database you may have to update documents preexisting in the database to reflect these changes**. To add a record, simply add the record name as a variable assignment and then assign the variable to a [MongoEngine field object](https://docs.mongoengine.org/guide/defining-documents.html#fields).