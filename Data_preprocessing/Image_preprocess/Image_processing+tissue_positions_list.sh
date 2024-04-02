### 1. need to run in the spatial environment
### 2. need to put spatial.txt file in the same directory


in_image="Spnb70_ROI.png"  # png image file
img_resolution_smaller=1300 # change this based on image resolutions, use smaller resolution


###
block_size_threshold=35   # default value is 35. block size for adaptive thresholding. This value may need adjusting depending on the specific characteristics of the image
C_normalized_in_tissue=0.1   # default value is 0.1. normalized value to determine the in tissue pixel (0-1)


python Image_processing+tissue_positions_list.py \
--img_path $in_image \
--block_size $block_size_threshold --C_normalized $C_normalized_in_tissue --img_res $img_resolution_smaller

cp $in_image ./spatial/tissue_lowres_image.png