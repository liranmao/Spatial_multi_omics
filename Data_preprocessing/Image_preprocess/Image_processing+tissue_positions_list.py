#!/usr/bin/env python
# coding: utf-8

# ## 1. Image Processing to generate binary image
import numpy as np
import matplotlib.pyplot as plt
from skimage.filters import threshold_local

from PIL import Image
import cv2
import pandas as pd
import imutils

import argparse

import math
import os

ap = argparse.ArgumentParser()
ap.add_argument("-i1", "--img_path", required=True, help="input image")
ap.add_argument("-i2", "--block_size", required=True, type=int, help="block size for thresholding")
ap.add_argument("-i3", "--C_normalized", required=True, type=float, help="normalized value to determine the in tissue pixel (0-1)")
ap.add_argument("-i4", "--img_res", required=True, type=int, help="image resolution (the smaller one)")

args = vars(ap.parse_args())


# Load image
img_path = args["img_path"]
img = Image.open(img_path)

# Convert image to grayscale if it's not
if len(np.array(img).shape) > 2:
    img = img.convert("L")

# Convert image to numpy array
img_array = np.array(img)

# Apply adaptive thresholding
block_size = args["block_size"] # default value is 35  # This value may need adjusting depending on the specific characteristics of the image
adaptive_thresh = threshold_local(img_array, block_size, offset=10)
binary_adaptive = img_array > adaptive_thresh

# Invert the image colors
binary_adaptive_inverted = np.invert(binary_adaptive)

# Convert numpy array back to PIL image and save
binary_adaptive_inverted_img = Image.fromarray(binary_adaptive_inverted.astype("uint8") * 255)

if not os.path.isfile('ROI_binary.png'):
    binary_adaptive_inverted_img.save("./ROI_binary.png")
else:
    print('ROI_binary.png exists. Will process the existing one.')


# ## 2. Generate_tissue_position_list

# Read image
I = cv2.imread('./ROI_binary.png', 0)

# Convert image to binary
_, BW = cv2.threshold(I, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)


pixel = 50
pixel_count = 2*pixel-1
numRows, numCols = BW.shape
pixel_w = numCols/pixel_count
pixel_h = numRows/pixel_count
str_data = ""

for i in range(1, 51):
    x = round(2*(i-1)*pixel_h + 1) - 1
    for j in sorted([i for i in range(1, 51)], reverse=True):
        y = round(2*(50-j)*pixel_w + 1) - 1
        img_pixel = BW[x:round(x+pixel_w-1), y:round(y+pixel_h-1)]
        C = np.sum(img_pixel)
        C_normalized = C / (255.0 * pixel_w * pixel_h) ### new line
        if C_normalized > args["C_normalized"]: # default value is 0.1 or use C > 0 if not use normalized intensity for thresholding
            str_data += f',{j}x{i}'


# Split the string by comma to get a list
list_values = str_data.split(',')

# Remove the first value and any empty strings from the list
list_values = [value for value in list_values[1:] if value]

# Convert the list into a DataFrame
position = pd.DataFrame(list_values, columns=['Values'])

# Add a new column 'in_tisue'
position['in_tisue'] = 1

# Set the row index to the values of the first column
position.index = position['Values']

# Read the spatial_barcodes.txt file
BC = pd.read_csv('./spatial_barcodes.txt', sep='\t', header=None)

# Create a new column 'position' by concatenating the values of the second and third columns
BC['position'] = BC[1].astype(str) + 'x' + BC[2].astype(str)

# Set the row index to the values of the 'position' column
BC.index = BC['position']

# Merge the two dataframes
df_in_tissue = pd.merge(position, BC, left_index=True, right_index=True, how='outer')

# Select only the required columns
df_in_tissue = df_in_tissue[[0, 'in_tisue', 2, 1]]


# Replace NaN values with 0
df_in_tissue = df_in_tissue.fillna(0)


# Subtract 1 from the values of the third and fourth columns
df_in_tissue[1] = df_in_tissue[1] - 1
df_in_tissue[2] = df_in_tissue[2] - 1





pixel_smaller = args["img_res"] # change this based on image resolutions
pixel_ratio = pixel_smaller / 198 


df_in_tissue[3] = (4*df_in_tissue[2]+1)*pixel_ratio
df_in_tissue[4] = (4*(49-df_in_tissue[1])+1)*pixel_ratio


df_in_tissue[3] = np.ceil(df_in_tissue[3])
df_in_tissue[4] = np.ceil(df_in_tissue[4])


dir_name = "./spatial/"

# Create the directory if it does not exist
if not os.path.exists(dir_name):
    os.makedirs(dir_name)

df_in_tissue.to_csv(dir_name+"./tissue_positions_list.csv", index=False, header=False)


# ## 3. Check_spots_alignment

im = cv2.imread(img_path)
im_res = im


spots = pd.read_csv(dir_name+"tissue_positions_list.csv", header=None)
spots = spots.set_axis(['barcode', 'in_tissue', 'array_row', 'array_column', 'pxl_col_in_fullres', 'pxl_row_in_fullres'], axis=1)


spot_width = pixel_smaller / 198
spot_height = pixel_smaller / 198


for i, spot in spots.iterrows():
    startX = int(round(spot['pxl_row_in_fullres'] - spot_width))
    startY = int(round(spot['pxl_col_in_fullres'] - spot_height))
    width = int(round(spot_width)*2)
    height = int(round(spot_height)*2)
    
    if(spot['in_tissue'] == 1):
        cv2.rectangle(im_res, (startX, startY), (startX+width, startY+height), (0, 255, 0), 1)
    #else:
        #cv2.rectangle(im_res, (startX, startY), (startX+width, startY+height), (0, 0, 255), 1)
        
cv2.imwrite('tissue_lowres_image_align_check.png', im_res) 


# ### generate json
# generate scalefactors_json.json
import json
scalefactors = {"spot_diameter_fullres": spot_width, 
                "tissue_hires_scalef": 1.0, 
                "fiducial_diameter_fullres": spot_width, 
                "tissue_lowres_scalef": 1.0}

with open(dir_name+'scalefactors_json.json', 'w') as outfile:
    json.dump(scalefactors, outfile)

print("Done!")