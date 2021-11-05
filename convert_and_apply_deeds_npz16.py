import numpy as np
import nibabel as nib
import struct
import scipy.ndimage
from scipy.ndimage.interpolation import zoom as zoom
from scipy.ndimage.interpolation import map_coordinates

import torch
import torch.nn as nn
import torch.nn.functional as F



import argparse

def transform(moving,disp_field):
    
        #disp_field = self.numpy_loader.load_image(disp_field_path).astype('float32')
        #upsample
        #disp_field = np.array([zoom(disp_field[i], 2, order=2) for i in range(3)])
        
        D, H, W = fixed.shape
        identity = np.meshgrid(np.arange(D), np.arange(H), np.arange(W), indexing='ij')
        moving_warped = map_coordinates(moving, identity + disp_field, order=0)


def main():

    parser = argparse.ArgumentParser()
    #inputdatagroup = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("--inputdat", dest="inputdat", help="input deeds displacement from (.dat)", default=None, required=True)
    parser.add_argument("--inputmatrix", dest="inputmatrix", help="input linear matrix file (.txt)", default=None, required=True)
    parser.add_argument("--inputseg", dest="inputseg", help="input segmentation (.nii.gz)", default=None, required=True)

    parser.add_argument("--outputseg", dest="outputseg",  help="output segmentation (.nii.gz)", default=None, required=True)


    options = parser.parse_args()
    d_options = vars(options)
   
    with open(d_options['inputdat'], 'rb') as content_file:
        content = content_file.read()
    H = 192
    W = 192
    D = 208
    
    grid_x = torch.arange(H).float().view(-1,1,1).repeat(1,W,D)
    grid_y = torch.arange(W).float().view(1,-1,1).repeat(H,1,D)
    grid_z = torch.arange(D).float().view(1,1,-1).repeat(H,W,1)
    transform_grid = torch.stack((grid_x,grid_y,grid_z),3).unsqueeze(0)

    grid_space = int((torch.pow(torch.Tensor([H*W*D])/(len(content)/12),0.334)))
    disp_field = torch.from_numpy(np.array(struct.unpack('f'*(len(content)//4),content))).reshape(1,3,D//grid_space,W//grid_space,H//grid_space).permute(0,1,4,3,2).float()
    disp_field = F.interpolate(disp_field,size=(H,W,D),mode='trilinear',align_corners=None).permute(0,2,3,4,1)[:,:,:,:,torch.Tensor([2,0,1]).long()].flip(4)
    #transform_grid+ #will be added later
    x = disp_field[0,:,:,:,0].numpy()
    y = disp_field[0,:,:,:,1].numpy()
    z = disp_field[0,:,:,:,2].numpy()
    
    
    D, H, W = fixed.shape
    identity = np.meshgrid(np.arange(D), np.arange(H), np.arange(W), indexing='ij')
    moving_warped = map_coordinates(moving, identity + disp_field, order=0)

    #x1 = zoom(x,1/2,order=2).astype('float16')
    #y1 = zoom(y,1/2,order=2).astype('float16')
    #z1 = zoom(z,1/2,order=2).astype('float16')

    

    #np.savez_compressed(d_options['outputnpz'],np.stack((x1,y1,z1),0))


if __name__ == '__main__':
    main()
