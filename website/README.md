# CIPHER Website Usage Guide

Welcome to the CIPHER database website usage guide, where you will learn how to setup, use, and modify the web application which is used to access the CIPHER database. Please note that all modification and testing described in this document is done on a locally hosted development server. Before deploying the original or modified application to a production server, please check to ensure that any development settings have been deactivated/changed to production and the website backend is properly secured and appropriate for the production server you are using. 

## Prerequisites and Setup

The web application frontend it written in JavaScript and the backend is written in Python using Flask. In order to run the website locally for testing, you must have a Conda environment active with the following packages installed:

```(yaml)
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

These packages are listed in the local `requirements.txt` file for convenient installation via Pip or Anaconda. To install the packages please follow the instructions below (note this assumes that the cwd is `cipher_db/website`):

### Pip

`pip install -r requirements.txt`

### Anaconda

Create a new conda environment if one for does not already exist for CIPHER:

`conda create -n python==3.7 cipher_website pip`

Activate the new enviroment:

`source activate cipher_website`

Install the requirements file via pip:

`pip install -r requirements.txt`

## Running the Web Application Locally

With the conda environment active, navigate to the `website` folder and run the `app.py` python file. From here, you should see the following text in the terminal:

```
* Serving Flask app 'app' (lazy loading)
* Environment: production
  WARNING: This is a development server. Do not use it in a production deployment.
  Use a production WSGI server instead.
* Debug mode: off
* Running on http://127.0.0.1:{PORT_NUMBER}/ (Press CTRL+C to quit)
```
Either `ctrl+click` on the link or type address with the corresponding port number into a web browser to view the web application.

**Note:** Flask will attempt to find an available port on the localhost automatically. The port number which the Flask app runs on can be set by going to line 170-171 of `app.py` and changing the following line of code with `{PORT_NUMBER}` set to the desired port. This should be changed prior to deploying on a development server.

```(python)
if __name__ == "__main__":
    app.run(host='localhost', port={PORT_NUMBER})
```

**Note:** If `app.py` is run on a remote server for local debugging, the port used will need to be forwarded from the remote server to the local host in order to view the webpage locally.

## Website Guide

### Homepage

The CipherDB homepage has 3 main sections which are used to provide introductory information and the current top ranked non-addictive analgesic compounds. The page looks as follows: 

![CipherDB Web Application Home Page](/cipher_db/website/README_images/Cipher_Webapp_Home.png)

The 'About CipherDB' section has general information about the project goals and how they are accomplished. The upper figure on the right hand side is the main overview figure for the project. The chart on the lower right hand side is an interactive chart which shows the top ranked non-addictive analgesics in the database at the current time. This is updated via a scheduled trigger to evolve over time as new non-addictive compounds are added to the database.

### Search Page

The CipherDB search page is used to search the database for compounds of interest. Compounds can be searched by a number of identifiers including name, SMILES, and InChI Key. A feature currently under development is to allow for compound searching based on chemical properties or biosignature results. Once the compound search information has been inputted, press the search icon on the right hand side and a summary of all compounds which match the description term will appear. This is shown below:

![CipherDB Web Application Search Page](/cipher_db/website/README_images/Cipher_Webapp_Search.png)

After searching, click on one of the four "snapshot" options ('Properties','Biosignature','Assays', and 'Synthesis') to view a simplified overview of that category for the given compound. Below is an example of the 'Properties' and 'Assay 'overview for Naloxone:

Properties

![CipherDB Web Application Properties Snapshot](/cipher_db/website/README_images/Cipher_Webapp_Properties_Snapshot.png)

Assays

![CipherDB Web Application Assay Snapshot](/cipher_db/website/README_images/Cipher_Webapp_Assay_Snapshot.png)

### Summary Page

Clicking on a compound on the search page will navigate to the summary page for that compound. On this page, the same four "snapshot" options will be available ('Properties','Biosignature','Assays', and 'Synthesis'), but clicking on them will provided a more detailed overview of the selected category. An example of this summary for Naloxone 'Properties' and 'Assays' is shown below:

Properties

![CipherDB Web Application Properties Summary](/cipher_db/website/README_images/Cipher_Webapp_Properties_Summary.png)

Assays

![CipherDB Web Application Assays Summary](/cipher_db/website/README_images/Cipher_Webapp_Assays_Summary.png)

### Add Page

The add page is designed to provide a GUI for adding compounds of interest to the database. This makes it simple for individuals not familiar with CLI-based methods to easily contribute novel compounds of interest to the database. An example of the add page is shown below:

**Insert Image Here**

The add page functions by inputting a SMILES string and designated name for the compound (optional) and from here the series of database triggers will mine (if applicable) and calculate chemical properties, determine CANDO binding information, and attempt to find a retrosynthetic route for the compound.

### Contact Page

This page contains contact information for the authors of this project.

### REST API

Fill this out later ...

## Website Files

### `app.py` - The Backend Flask Application

This python file contains the main backend Flask code for the web application. Each python method corresponds to a specific URL and based on the request type, performs the correspond functions. Please note that only GET and POST requests are currently supported by the backend. The following methods are currently implemented:
| Method    | Route | Description |
| ----------| ----- | ----------- |
| `index()` | `/`   | GET request only. Renders the `index.html` page upon GET request |
| `add()`   | `/add`| Renders the `add.html` page upon GET request. Access the inputted compound name and SMILES and add to the database (with trigger support) upon POST request. |
| `top()`   | `/top`| POST request only. Used for access to the top ranked compounds based on CANDO biosignature. |
| `search()`| `/search` | Renders search page template upon GET request. Upon POST request, sends the search information to the website backend, queries the database, and returns compounds matching the search term.
| `info()`  | `/info`| POST request only. Returns compound snapshot information for quick summary as shown on the search page|
| `summarty()` | `/summary/<inchikey>` | GET request only. Returns compound summary information from the database for the specified compound as shown on the compound summary page.
| `contact()` | `/contact` | Fill out later based on implementation ... |
| `rest_*()`  | `rest/...` | REST methods. See the REST API section for more details |

### `engine.py` - The Backend Database Engine

This python file contains the scripts which facilitate interactions between the website backend and the database. An overview of the methods are given below:

| Method | Arguments | Description |
| ------ | --------- | ----------- |
| `return_compound()`| `identifier` | Queries the `Compounds` collection of the database for the given identifier and returns JSON formatted information on the selected compounds.|
| `return_properties()` | `inchikey` | Queries the `Properties` collection of the database for the compound with the given identifier and returns JSON formatted information on the compounds properties as mined from PubChem and calculated via RDKit.|
| `return_biosignature()` | `inchikey`| Queries the `Assay` collection of the database for literature derived MOA information for the compound. Queries the `CANDO` collection of the database for biosignautre information for the selected compound. Returns the information in python dictionary format.|
| `return_assay()` | `inchikey`| Queries the `Assay` collection of the database for literature derived assays from PubChem, Drugbank, CHEMBL, and the IUPHAR database, and the DESI MS instrumented used in the ASPIRE project on the specific compound. Returns the Assay information in JSON format.|
| `return_askcos_pathway()` | `inchikey` | Queries the `Retorsyntheisis` collection of the database for retrosynthetic pathway information on the specified compound. Returns a list of images of the retrosynthetic pathways.|
| `return_compound_image()` | `inchikey` | Queries the `Compounds` collection of the database for the compound image as calculated using RDKit. Returns the image of the specified compound.|
| `return_desired_dynamic_biosignature()` | `None` | Returns the desired dynamic biosignature as specified in the aims of the ASPIRE project in python dictionary format.|
| `return_biosig_knn()` | `None` | Queries the `CANDO` collection of the database for the k nearest neighboring biosignatures of the desired dynamic biosignature. Returns the k nearest biosignatures as a list of python dictionaries.|


### `passenger_wsgi.py` - Runnable Web Application

This file contains a runnable version of the Flask application for web server deployment. The code is shown below:

```(python)
from app import app as application
```

When setting up a python web application which needs to call the runnable Flask application, make sure to specify the file is `passenger_wsgi.py` and the application is called `application`.

### `templates/new/*.html` - Frontend HTML Webpage Templates

This folder contains the HTML templates for the web application.

- `home.html`: Homepage as shown in the website guide
- `search.html`: Search page as shown in the website guide
- `summary.html`: Summary page as shown in the website guide
- `add.html`: Add page as shown in the website guide

### `static/assets/js/*.js` - Javascript Files for Website Frontend

This folder contains the Javascpript files for the web application.

- `home.js`: Retrieves the desired biosignature and top 5 most similar biosignatures from the backend. These are then rendered in an interactive heatmap-like figure on the home page. An overview of the methods are given below:

| Method | Arguments | Description |
| ------ | --------- | ----------- |
| `getList()`| `None` | Retrieves biosignatures from the backend and calls `renderResults()` |
| `renderBiosigs()` | `desired`,   `top` | Uses the ApexCharts library to render the interactive biosignature figure with the `desired` biosignature and the `top` 5 similar biosignatures. |
| `renderResults()` | `results`| Formats backend results to allow for proper rendering. Calls `renderBiosigs()`  |

- `search.js`: Manages communication between the frontend and backend on the search page. Retrieves results from the backend and renders them for the user. An overview of the methods are given below:

| Method | Arguments | Description |
| ------ | --------- | ----------- |
| `capitalizeFirstLetter()` | `string` | Capitalizes the first letter of the argument and returns the result. Used when rendering compound names retrieved from the backend. |
| `search()` | `None` | Performs a search request with the term currently entered in the search bar. Calls `renderResults()` with the response provided by the backend. |
| `load()` | `query` | Renders a loading animation while waiting for `search()` to finish. |
| `displayToggle()`| `cmpdName`, `toggle` | Manages which content is displayed in the snapshot view under a given search result corresponding to `cmpdName`. | 
| `renderBiosig()` | `cmpdName`, `inchikey`, `cmpdBiosig`, `desired` | Renders the biosignature of a given search result corresponding to `inchikey`. Renders the `desired` biosignature below this for comparison. | 
|  `renderResults()` | `results` | Manages the rendering of all search results and snapshot views. | 

- `summary.js`: Retrieves summary information from the backend for a given compound and renders this on the `summary.html` page. An overview of the methods are given below:

| Method | Arguments | Description |
| ------ | --------- | ----------- |
| `capitalizeFirstLetter()` | `string` | Capitalizes the first letter of the argument and returns the result. Used when rendering compound names retrieved from the backend. |
| `search()` | `query` | Retrieves compound summary information corresponding to the inchikey provided in the url. Calls `renderResults()` with the response provided by the backend. |
| `renderBiosig()` | `cmpdName`, `inchikey`, `cmpdBiosig`, `desired` | Renders the biosignature of a given search result corresponding to `inchikey`. Renders the `desired` biosignature below this for comparison. | 
|  `renderResults()` | `results` | Manages the rendering of current summary information and snapshot views. | 

- `add.js`: Manages communication between the frontend and backend on the `add.html` page. An overview of the methods are given below:

| Method | Arguments | Description |
| ------ | --------- | ----------- |
| `add()` | `query` | Submits the SMILES string present in the text input on `add.html` to the backend for processing. |

- `bs-init.js`: Manages tooltips and responsive pages.  

## Modifications FAQ

Modifications for the web application are relatively simple provided the modifier has basic to intermediate web development experience. Since the modification scope for this project can be very broad, this will provide a general overview of areas where modifications may need to be made if new fields/collections are added to the database and/or new modules are added. No website frontend modifications are covered, all modifications are either for the website backend or database.

#### Question
"I added a new type of identifier I would like be able to search by"

#### Answer
Go to the `return_compounds()` method of in `engine.py` and include an `elif` statement with the following text (replacing `<identifier_name>`):

```(python)
elif Compounds.objects(<identifier_name>=identifier).count() > 0:
  return json.loads(Compounds.objects(<identifier_name>=identifier).to_json())
```

#### Question

"I want to change the proteins queried in the biosignature"

#### Answer

Go to the `return_biosignature()` method in `engine.py` and modify the `BIOSIG` and `BMIDS` constants to include the name and 6 character biomolecule ID in the appropriate position of the list and dictionary.

#### Question

"I want to include a new data source for MOA information for biosignatures which I have added to the database `Assay` collection"

#### Answer

Go to the `return_biosignature()` method in `engine.py` and add a constant `<DATASOURCE_NAME>` with the name of your data source. Then add a series of `if` statements in the `for` loop formatted in the following manner:

```(python)

```
