import numpy as np
import scipy.sparse as sp
import torch
from os import listdir
from os.path import isfile, join
import scipy.io as sio
from tqdm import tqdm
from utils.mis_functions.graph_sample import GraphSample

data_name = 'mindboggle'
data_path = './output/' + data_name +'_forpy/'
save_path = './output/pyt_files_' + data_name 

if not os.path.exists(save_path):
  os.makedirs(save_path)

files = [f for f in listdir(data_path) if isfile(join(data_path, f))]

for cnt in tqdm(range(len(files))):

    temp = sio.loadmat(data_path + '/' + files[cnt])
    x = torch.cat((torch.FloatTensor(temp['X'][:, 0:3]),
                   torch.FloatTensor(temp['C']),
                   torch.FloatTensor(temp['T'])), 1)
    e1, e2, e3 = sp.find(temp['A'])
    edge_idx = torch.cat((torch.LongTensor(e2).unsqueeze(0), torch.LongTensor(e1).unsqueeze(0)),
                         0)  # index x,y of the sparse adj matrix
    edge_wht = torch.FloatTensor(e3).unsqueeze(0)  # the weights of the edges. Often used to construct the sparse adj
    gt = torch.LongTensor(temp['GT'])  # This is manual parcel label for mindboggle;
    xyz = torch.FloatTensor(temp['EUC'])  # This is the xyz of the mesh node location in euclidean coordinates
    age = torch.FloatTensor(temp['AG'][0])  # This is the age of the subject
    sx = torch.FloatTensor(temp['SX'])  # This is the gender of the subject
    lab = torch.FloatTensor(temp['Y'])  # This is FreeSurfer label
    data = GraphSample(x=x, edge_idx=edge_idx, edge_wht=edge_wht, gt=gt, xyz=xyz,
                       age=age, sx=sx,
                       lab=lab)
    torch.save(data, save_path + files[cnt][:-4] + '.pt')
