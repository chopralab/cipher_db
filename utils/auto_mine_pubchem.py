import argparse
from io import StringIO
import json
import os
from typing import Dict, List, Optional, Union
import xml.etree.ElementTree as ET
import warnings

from PIL import Image
import pandas as pd
import requests

# Exception for invalid command line arguments
class ArgumentError(Exception):
    pass

# Exception for PubChem REST API error status codes
class StatusCodeError(Exception):
    def __init__(self, status_code):
        message = f'Unknown Status Code Received: {status_code}'

        if status_code == 400:
            message = 'Status Code 400 - PUGREST.BadRequest - Request is improperly formed (syntax error in the URL, POST body, etc.)'
        elif status_code == 404:
            message = 'Status Code 404 - PUGREST.NotFound - The input record was not found (e.g. invalid CID)'
        elif status_code == 405:
            message = 'Status Code 405 - PUGREST.NotAllowed - Request not allowed (such as invalid MIME type in the HTTP Accept header)'
        elif status_code == 504:
            message = 'Status Code 504 - PUGREST.Timeout - The request timed out, from server overload or too broad a request'
        elif status_code == 503:
            message = 'Status Code 503 - PUGREST.ServerBusy - Too many requests or server is busy, retry later'
        elif status_code == 501:
            message = 'Status Code 501 - PUGREST.Unimplemented - The requested operation has not (yet) been implemented by the server'
        elif status_code == 500:
            message = 'Status Code 500 - PUGREST.ServerError - Some problem on the server side (such as a database server down, etc.)'
 
        super().__init__(message)

# Set up the argument parser
parser = argparse.ArgumentParser()

# General setup
parser.add_argument('--infile', type=str, default='infile.txt', help='The name/path of the file containing input information to the REST API')
parser.add_argument('--outfile', type=str, default='outfile.txt', help='The name/path of the file that will contain the data output of the REST API')
parser.add_argument('--no-append', action='store_true', default=False, help='Flag for creating a seperate file for the output of each REST API call')

# If you just want to mine a list of urls
parser.add_argument('--url', action='store_true', default=False, help='Flag for URL mode; select this for running a file of PubChem REST API URL\'s directly')

# Override for REST API domain
parser.add_argument('--domain-override', choices=['substance','compound','assay'], help='Override the domain of access (Don\'t use this unless you are certian you need to change domains)')

# Pubchem REST API input
parser.add_argument('--in-cid', action='store_true', default=False, help='Input - Flag for compound ID')
parser.add_argument('--in-sid', action='store_true', default=False, help='Input - Flag for substance ID')
parser.add_argument('--in-name', action='store_true', default=False, help='Input - Flag for compound name')
parser.add_argument('--in-smiles', action='store_true', default=False, help='Input - Flag for compound SMILES string')
parser.add_argument('--in-formula', action='store_true', default=False, help='Input - Flag for compound forumla')
parser.add_argument('--in-inchikey', action='store_true', default=False, help='Input - Flag for compound InChi Key')
parser.add_argument('--in-xref', action='store_true', default=False, help='Input - Flag for cross reference identifiers, must be used with --in-xtype (i.e. serach by patent identification number)')
parser.add_argument('--in-xtype', choices=['RegistryID','RN','PubMedID','MMDBID','ProteinGI','NucleotideGI','TaxonomyID','MIMID','GeneID','ProbeID','PatentID'], help='Input - Cross reference type, must be used with --in-xref (i.e. \"PatentID\")')
parser.add_argument('--in-structure', choices=['substructure','superstructure','similarity','identity'], help='Input - Flag for special structure search (must be used with --in-smiles or --in-cid)')
parser.add_argument('--direct-input', action='store_true', default=False, help='Input - Argument for providing input URL segment directly (i.e. \"compound/name/glucose/\")')

# Pubchem REST API data
parser.add_argument('--data-record', action='store_true', default=False, help='Data - Flag for full PubChem record of the selected compound')
parser.add_argument('--data-image', action='store_true', default=False, help='Data - Flag for PNG image of the selected compound')
parser.add_argument('--data-property', type=str, help='Data - Argument for chemical property (or properties) of the selected compound (for multiple properties do \" MolecularWeight,MolecularFormulat,InChi,...\")')
parser.add_argument('--data-synonyms', action='store_true', default=False, help='Data - Flag for synonyms of the selected compound')
parser.add_argument('--data-xref', action='store_true', default=False, help='Data - Flag for cross reference identifiers, must be used with --in-xtype (i.e. serach by patent identification number)')
parser.add_argument('--data-xtype', choices=['RegistryID','RN','PubMedID','MMDBID','ProteinGI','NucleotideGI','TaxonomyID','MIMID','GeneID','ProbeID','PatentID'], help='Data - Cross reference type, must be used with --in-xref (i.e. \"PatentID\")')

# Pubchem BioAssay input/data
parser.add_argument('--assay', action='store_true', default=False, help='BioAssay - Flag for use of PubChem BioAssay database')
parser.add_argument('--aid', action='store_true', help='BioAssay Input - Flag for assay ID')
parser.add_argument('--atarget', type=str, help='BioAssay Data - Argument for assay data (i.e. \"ProteinGI,ProteinName,GeneID,...\")')
parser.add_argument('--adescription', action='store_true', default=False, help='BioAssay Data - Flag for getting assay description')
parser.add_argument('--asummary', action='store_true', default=False, help='BioAssay Data - Flag for getting assay summary')

# Pubchem REST API output
parser.add_argument('--outtype', choices=['XML','JSON','JSONP','ASNB','ASNT','SDF','CSV','PNG','TXT'], default='TXT', help='Output - Format of the output file')
parser.add_argument('--callback', type=str, help='Output - Used with JSONP output type')
parser.add_argument('--include-label', action='store_true', default=False, help='Output - Flag for including input compound in output file, included as @index - compound in text files above the information. Only works for text and CSV file output.')

# Parser the arguments
args = parser.parse_args()

class PubChem_Miner():

    def __init__(self):
        pass

    @staticmethod
    def get_compound_info(
        input_type: str, input: str, data_type: str,
        data: Union[List, str], output_type: str = 'TXT'
    ) -> Optional[Union[str, pd.DataFrame, Dict, ET.ElementTree]]:
        """
        get the specified data for a given compound

        Parameters
        ----------
        input_type : str
            the type of the input query
        input : str
            the compound for which to retrieve data for
        data_type : str
            the type of data to query for
        data : Union[List, str]
            the specific data to retrieve
        output_type : str, default='TXT'
            how the results fo the query should be returned

        Returns
        -------
        Optional[Union[str, pd.DataFrame, Dict, ET.ElementTree]]
            the output, if a valid query was made

        Raises
        ------
        StatusCodeError
            [description]
        """
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
        domain = "compound/"

        if input_type == 'smiles':
            input = 'smiles/' + input + '/'
        elif input_type == 'inchikey':
            input = 'inchikey/' + input + '/'
        elif input_type == 'name':
            input = 'name/' + input + '/'
        else:
            print("Input type not supported")
            print(type(input_type))
            print(input_type)
            print(type(input))
            print(input)
            return None

        if data_type is 'property':
            if isinstance(data, list):
                # dcopy = 'property/'
                data = 'property/' + ','.join(
                    ['canonicalSMILES' if d=='smiles' else d for d in data]
                ) + '/'
                # for elem in data:
                #     if elem == 'smiles':
                #         dcopy += 'canonicalSMILES' + ','
                #     else:
                #         dcopy += elem + ','
                # dcopy = dcopy[:-1]
                # data = dcopy + '/'
            else:
                if data == 'smiles':
                    data = 'property/' + 'canonicalSMILES' + '/'
                else:
                    data = 'property/' + data + '/'
        elif data_type == 'synonyms':
            data = 'synonyms' + '/'
        elif data_type == 'assay':
            if output_type not in ['XML','CSV','JSON']:
                print("Invalid output type for this query")
                return None
            else:
                data = 'assaysummary/' 
        else:
            print("Data type not supported")
            return None

        if output_type not in ['TXT','XML','CSV','JSON']:
            print('Output type not supported')
            return None

        url = url + domain + input + data + output_type
        response = requests.get(url)
        if response:
            print("Request on " + url + " succeded!")
            if output_type == 'TXT':
                return response.text
            if output_type == 'CSV':
                return pd.read_csv(StringIO(response.text))
            if output_type == 'XML':
                return ET.ElementTree(ET.fromstring(StringIO(response.text)))
            if output_type == 'JSON':
                return json.loads(StringIO(response.text))
            
            raise ValueError(f'Invalid output_type: "{output_type}"')
        else:
            try:
                raise StatusCodeError(response.status_code)
            except StatusCodeError:
                print("Request on " + url + " failed")
                return None

# Check to see if all argument values are satisfied
def check_args():
    num_warnings = 0
    if not os.path.isfile(args.infile):
       raise ArgumentError("Infile not found (type --help for more information)")

    if args.url is False:
        if args.outtype is 'JSONP' and args.callback is None:
            warnings.warn("Output type JSONP is used with no callback type specified, default \"callback\" will be used")
            num_warnings += 1
        if args.data_property is not None and len(args.data_property.split(",")) > 1 and args.outtype == 'TXT':
            num_warnings += 1
            warnings.warn("Cannot use output type \"TXT\" when requesting more than one chemical property, defaulted to CSV")
            args.outtype = 'CSV'   
    if args.url is False and args.assay is False:
        if (args.in_structure is not None and args.in_smiles is False) and (args.in_structure != None and args.in_cid is False):
            raise ArgumentError("Cannot use --in-structure without --in-smiles or --in-cid (type --help for more information)")
        if args.in_xref and args.in_xtype is None:
            raise ArgumentError("Cannot use --in_xref without specifying and xref type using --in_xtype (type --help for more information)")
        if args.in_xref is False and args.in_xtype is not None:
            warnings.warn("Argument --in_xtype will not be used due to not flagging --in_xref")
            num_warnings += 1
    
        num_input = 0
        if args.in_cid:
            num_input += 1
        if args.in_sid:
            num_input += 1
        if args.in_name:
            num_input += 1
        if args.in_smiles:
            num_input +=1
        if args.in_formula:
            num_input += 1
        if args.in_inchikey:
            num_input += 1
        if args.in_xref:
            num_input += 1
        if args.direct_input:
            num_input += 1
        
        if num_input == 0:
            raise ArgumentError("An input flag is needed to determine the type of input (type --help for more information)")
        if num_input > 1:
            raise ArgumentError("More than one input flag is not allowed (type --help for more information)")
          
        num_data = 0
        if args.data_record:
            num_data += 1
        if args.data_image:
            num_data += 1
        if args.data_property is not None:
            num_data += 1
        if args.data_synonyms:
            num_data += 1
        if args.data_xref:
            num_data += 1
        
        if num_data == 0:
            args.data_record = True
            warnings.warn("No data flag provided, default --data_record is selected (type --help for more information)")
            num_warnings += 1
        
        if args.data_xref and args.data_xtype is None:
            raise ArgumentError("Cannot use --data_xref without specifying and xref type using --data_xtype (type --help for more information)")
        if args.data_xref is False and args.data_xtype is not None:
            warnings.warn("Argument --data_xtype will not be used due to not flagging --data_xref")
            num_warnings += 1
    if args.assay:
            if args.aid is False and args.in_cid is False and args.in_sid is False:
                raise ArgumentError("If assay is selected either --cid, --sid, or --aid must be provided (type --help for more information)")

            if args.atarget is None and args.adescription is False and args.asummary is False:
                args.asummary = True
                warnings.warn("No bioassay data is selected, default summary will be selected")
                num_warnings += 1
    
    if num_warnings > 0:
        print("Arguments parsed with " + str(num_warnings) + " warning(s)")
    else:
        print("Arguments parsed with no errors or warnings")

# Construct a PubChem REST API URL from given information
def construct_url(line):
    # Initalize the URL components
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
    domain = ""
    input = ""
    data = ""
    output = ""

    # Assign the domain and input
    if args.in_structure is not None:
        input += args.in_structure + "/"
    elif args.in_cid:
        domain = "compound/"
        input += "cid/" + line + "/"
    elif args.in_sid:
        domain = "substance/"
        input += "sid/" + line + "/"
    elif args.in_name:
        domain = "compound/"
        input += "name/" + line + "/"
    elif args.in_formula:
        domain = "compound/"
        input += "formula/" + line + "/"
    elif args.in_inchikey:
        domain = "compound/"
        input += "inchikey/" + line + "/"
    elif args.in_xref:
        domain = "substance/"
        input += "xref/" + args.in_xtype + "/" + line + "/"
    elif args.direct_input:
        domain = "substance/"
        input = line
    
    # Override the domain if requested
    if args.domain_override is not None:
        domain = args.domain_override + "/"

    # Assign the data
    if args.data_record:
        data = "record/"
    elif args.data_image:
        args.outtype = "PNG"
    elif args.data_property is not None:
        data = "property/" + args.data_property + "/"
    elif args.data_synonyms:
        data = "synonyms/"
    elif args.data_xref:
        data = "xrefs/" + args.data_xtype + "/"
    
    # Assignment for assay data
    if args.aid:
        domain = "assay/"
        input = "aid/" + line + "/"
    if args.atarget is not None:
        data = "targets/" + args.atarget + "/"
    if args.adescription:
        data = "description/"
    if args.asummary:
        data = "summary/"

    # Assign the output
    output = args.outtype
    if args.outtype == "JSONP" and args.callback is not None:
        output += "?callback=" + args.callback

    url = url + domain + input + data + output
    return url

# Call a PubChem REST API URL and receive a respone
def pugrest_request(url):
    response = requests.get(url)
    if response:
        print("Request on " + url + " succeded!")
        return response
    else:
        try:
            raise StatusCodeError(response.status_code)
        except StatusCodeError:
            print("Request on " + url + " failed")
            return response

# Reconstruct the output url for no append mode
def construct_out_url(line, url):
    url_parts = url.split("/")
    out_url = ""
    for i in range(len(url_parts)):
        if i == 0:
            out_url += url_parts[i]
        elif i == len(url_parts)-1:
            out_url += "/" + line + "_" + url_parts[i]
        else:
            out_url += "/" + url_parts[i]
    return out_url

# Save the PubChem REST API output in the format specified by the command line arguments
def save_output(outtype, response, line):
    out_url = construct_out_url(line,args.outfile)
    if outtype == "PNG":
        with open(out_url,'wb') as outfile:
            outfile.write(response.content)
    elif outtype == "TXT":
        if not args.no_append:
            outfile = open(args.outfile,'a')
        else:
            outfile = open(out_url,'w')
        
        if args.include_label:
            outfile.write("@index - " + line + "\n")
            outfile.write(response.text)
        else:
            outfile.write(response.text)
        outfile.close()
    elif outtype == "CSV":
        df = pd.read_csv(StringIO(response.text))
        zero_row = [line] * len(df)
        zero_row = pd.DataFrame(zero_row, index=False)
        if not args.no_append:
            if args.include_label:
                df = zero_row.join(df)
                df.to_csv(args.outfile, mode='a',index=False)
            else:
                df.to_csv(args.outfile, mode='a',index=False)
        else:
            if args.include_label:
                df = zero_row.join(df)
                df.to_csv(out_url, index=False)
            else:
                df.to_csv(out_url, index=False)
    elif outtype == "JSON":
        json_string = json.dumps(json.load(response.json,strict=False))
        if not args.no_append:
            with open(args.outfile,'a') as outfile:
                json.dump(json_string,outfile)
            outfile.close()
        else:
            with open(out_url,'w') as outfile:
                json.dump(json_string,outfile)
            outfile.close()
    elif outtype == "XML":
        if not args.no_append:
            with open(args.outfile,'ab') as outfile:
                outfile.write(ET.tostring(response.text))
            outfile.close()
        else:
            with open(out_url,'wb') as outfile:
                outfile.write(ET.tostring(response.text))
            outfile.close()
    else:
        if not args.no_append:
            outfile = open(args.outfile,'a')
        else:
            outfile = open(out_url,'w')
        outfile.write(response.text)
        outfile.close()
    
# Call PubChem REST API URLs for all lines in the file
def process_all():
    with open(args.infile,'r') as infile:
        for line in infile:
            line = line.replace("\n","")
            if args.url:
                response = pugrest_request(line)
                if response:
                    outtype = line.split("/")[-1]
                    save_output(outtype, response, line)
            else:
                url = construct_url(line)
                response = pugrest_request(url)
                if response:
                    save_output(args.outtype, response, line)
    infile.close()
                
# Main method           
def main():
    check_args()
    process_all()

# Run the main method
if __name__ == "__main__":
    main()   