#!/path/to/conda/env/bin/python3
import numpy as np
import argparse
import os
from rdkit import Chem
import pandas as pd
from sklearn.preprocessing import minmax_scale
from keras.layers import Input, Dense
from keras.models import Model, load_model, Sequential
from keras import initializers
import tensorflow as tf
from CVAE.model import CVAE
from multiprocessing import Process
from sklearn.metrics import mean_squared_error
from keras import backend as K
from google_drive_downloader import GoogleDriveDownloader as gdd

parser = argparse.ArgumentParser()
parser.add_argument('-t','--template', help="Template objective signature path.",required=True)
parser.add_argument('-d','--designs',help="Number of designs to be generated with de novo signature.",required=True)
parser.add_argument('-o','--output',help="Filename to output designed SMILES.",required=True)

parser.add_argument('--batch_size', help=argparse.SUPPRESS, type=int, default=128)
parser.add_argument('--latent_size', help=argparse.SUPPRESS, type=int, default=200)
parser.add_argument('--unit_size', help=argparse.SUPPRESS, type=int, default=512)
parser.add_argument('--n_rnn_layer', help=argparse.SUPPRESS, type=int, default=3)
parser.add_argument('--mean', help=argparse.SUPPRESS, type=float, default=0.0)
parser.add_argument('--stddev', help=argparse.SUPPRESS, type=float, default=1.0)
parser.add_argument('--num_prop', help=argparse.SUPPRESS, type=int, default=200)
parser.add_argument('--lr', help=argparse.SUPPRESS, type=float, default=0.0001)

args = parser.parse_args()

def get_models():
   if os.path.exists("./models") == False:
      print("Downloading models. Do not quit the application until finished.")
      gdd.download_file_from_google_drive(file_id='18f5pD12l8ttsSKNj0dVeFsekkFPcsMri',dest_path='./models.zip',unzip=True)
      print("Done.")
   else:
      None

def read_sig(path):
   file = open(path,'r')
   vals = file.readlines()
   vals = [float(v.split('\t')[1]) for v in vals]
   return vals

def root_mean_squared_error(y_true, y_pred):
        return K.sqrt(K.mean(K.square(y_pred - y_true)))

def extract_layers(main_model, start, end):
  new_model = Sequential()
  for i in range(start, end + 1):
    curr_layer = main_model.get_layer(index=i)
    new_model.add(curr_layer)
  return new_model

def encode(vals,ae_path):
    restore = load_model(ae_path,custom_objects={'root_mean_squared_error':root_mean_squared_error})
    enc = extract_layers(restore,0,10)
    print(enc.summary())
    sig = pd.DataFrame([vals])
    pred = enc.predict(sig,batch_size=1)[0].tolist()
    return pred

def get_smiles(tokenized, char):
   return "".join(map(lambda x: list(char)[x], tokenized.astype(int))).strip()

def design(encoded_sig,model_path):
    tf.reset_default_graph()
    char = ('c', 'C', '(', ')', '1', 'O', '2', '=', 'N', 'n', '3', 'S', '-', '[', ']', '/', 'l', 'F', 's', '4', '+', 'o', 'H', 'B', 'r', '#', '\\', 'I', '5', '@', 'P', '6', '7', 'E', 'X')
    vocab = {'c': 0, 'C': 1, '(': 2, ')': 3, '1': 4, 'O': 5, '2': 6, '=': 7, 'N': 8, 'n': 9, '3': 10, 'S': 11, '-': 12, '[': 13, ']': 14, '/': 15, 'l': 16, 'F': 17, 's': 18, '4': 19, '+': 20, 'o': 21, 'H': 22, 'B': 23, 'r': 24, '#': 25, '\\': 26, 'I': 27, '5': 28, '@': 29, 'P': 30, '6': 31, '7': 32, 'E': 33, 'X': 34}
    model = CVAE(len(char),args)
    model.restore(model_path)
    repeated_sig = np.array([encoded_sig for _ in range(128)])
    init_char = np.array([np.array(list(map(vocab.get, 'X'))) for _ in range(128)])
    smiles = []
    for _ in range(100):
        batch = model.sample(np.random.normal(0.0, 1.0, (128, 200)), repeated_sig, init_char, 120)
        smiles += [get_smiles(batch[i], char) for i in range(len(batch))]
    smiles = list(set([s.split('E')[0] for s in smiles]))
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    mols = [mol for mol in mols if mol is not None]
    smiles = [Chem.MolToSmiles(mol) for mol in mols]
    try:
       smiles = smiles[:5000]
    except:
       smiles = smiles
    return smiles

if __name__ == "__main__":
   get_models()
   template = args.template
   template_sig = read_sig(template)
   designs = []
   full = False
   while not full:
      for s in design(encode(template_sig,'models/ae.h5'),'models/model_30.ckpt-30'):
         designs.append(s)
      if len(designs) >= int(args.designs):
         designs = designs[:int(args.designs)]
         full = True
   designs = [d+'\n' for d in designs]

   out = open(args.output,'w')
   out.writelines(designs)
   out.close()