# Assays

## Background

In this project, the Assay collection of the database holds documents which pertain to a range of biological assays which were mined from various online databases as well as run in lab using the DESI-MS instrumentation. We have mined numerous bioassays relating to current opioid and opioid like compounds which have specific interactions with proteins listed in the target biosignature. The assays have all been mined from PubChem, but come from external databases such as CHEMBl, Drugbank, and the IUPHAR database. Additionally, this collection is designed to receive assay data from the instrument in an automated manner via the REST API.

The collection is designed using a flexible document structure to accommodate for the wide range of fields which may be encountered in various assays. In technical terms, this means that the MongoEngine code which handles database communication will not restrict each document to having specific fields (unlike in other modules). Thus, the document structure for this collection only enforces fields for identification purposes (CIPHER assay ID, InChIKey, SMILES, data source, etc.). 

## Document Structure

| Variable Name | Type | Description |
| ------------- | ---- | ---------- |
| `cipher_aid` | `str` | A unique identifier assigned to the assay |
| `inchikey` | `str` | The inchikey of the target compound of the assay | 
| `smiles` | `str` | The SMILES string of the target compound of the assay |
| `source` | `str` | The source of the assay |
| `receptor` | `str` | The target receptor of the assay |
| `{dynamic_field}` | `{dynamic_type}` | Any additional fields can be added due to the dynamic nature of the document |

## Standardization for Assays

Since there is no common structure to the assay document (besides identifiers), this can present challenges upon querying the collection for relevant information. One option to easily avoid this issue is to create an internal standard for field names for all assays which are added to the database, however this is not always feasible. Another potential option to help solve this issue is to assign grouping keywords or codes to felids (i.e. MOA_1, Mode_of_Action_1, MofA_1, etc.) and then when search for a field, you grab the keys for all fields of all documents and filter the search for fields with the keyword. A third option would be to use embedded documents to create objects for different assays.

## Modification

If you would like to change the document structure from a dynamic document to a regular document, simply add the required fields which you would like to enforce and then change the class definition to the following

```(python)
class Assays(me.Document)
```