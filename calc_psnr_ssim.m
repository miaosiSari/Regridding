function [p,s] = calc_psnr_ssim(path1, path2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  FUNCTION:  metrix_psnr
%%%
%%%  INPUTS:    reference_image     - original image data  path
%%%
%%%             query_image         - modified image data path to be compared with
%%%                                   original image
%%%
%%%  OUTPUTS:   psnr        - PSNR value
%%%
%%%  OUTPUTS:   ssim        - SSIM value
%%%
%%%  CHANGES:   NONE
%%%             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
img1 = im2uint8(imread(path1));
img2 = im2uint8(imread(path2));
p = metrix_psnr(img1, img2);
s = metrix_ssim(img1, img2);