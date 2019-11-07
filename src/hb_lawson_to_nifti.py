#!/usr/bin/env python

import nibabel as nib
import numpy as np
import os
import sys

filter_low_intensity = True

if len(sys.argv) < 3:
    sys.stderr.write('Usage: %s ref_nifti csv_filename [output_nifti_filename]\n' % os.path.basename(sys.argv[0]))
    sys.exit(-1)

fn_ref = sys.argv[1]
fn_csv = sys.argv[2]
if len(sys.argv) > 3:
    fn_out = sys.argv[3]
    if bn_out[-3:] == '.gz':
        bn_out = bn_out[:-3]
    if bn_out[-4:] == '.nii':
        bn_out = bn_out[:-4]
else:
    bn_out = fn_csv
    if bn_out[-4:] == '.csv':
        bn_out = bn_out[:-4]
    fn_out = bn_out + '.nii.gz'


def is_left(point, v0, v1):
    ax, ay = v0
    bx, by = v1
    cx, cy = point
    return ((ax-cx)*(by-cy) - (ay-cy)*(bx-cx)) > 0

def is_in(point, lst_pixel_ro):
    for i in range(len(lst_pixel_ro)-1):
        if not is_left(point, lst_pixel_ro[i], lst_pixel_ro[i+1]):
            return False
    return True

def add_nbd_to_set(point, set_pixel):
    x, y = round(point[0]+0.5), round(point[1]+0.5)
    set_pixel.add((x-0.5, y-0.5))
    
    return



def pixel_to_voxel(lst, lr='r'):
    check_left = (lr == 'r')
    voxel_y = int(lst[0])
    lst_pixel = [(lst[i], lst[i+1]) for i in range(1, len(lst), 2)]
    lst_pixel_ro = [lst_pixel[i] for i in [0, 3, 1, 4, 0]]
    
    if (0,0) in lst_pixel_ro:
        print '(0,0) in lst. skip: %s' % lst
        return []
    if lst[11] < 4:
        print '#pts < 4 in lst. skip: %s' % lst
        return []
    
    lst_xx = [pixel[0] for pixel in lst_pixel[:-1]]
    lst_yy = [pixel[1] for pixel in lst_pixel[:-1]]

    lst_xx_ro = [lst_xx[i] for i in [0, 3, 1, 4, 0]]
    lst_yy_ro = [lst_yy[i] for i in [0, 3, 1, 4, 0]]
    
    x0 = int(np.floor(min(lst_xx_ro)))
    xn = int(np.ceil(max(lst_xx_ro)))
    y0 = int(np.floor(min(lst_yy_ro)))
    yn = int(np.ceil(max(lst_yy_ro)))

    in_pixels = set()
    for x in np.arange(x0-1, xn+2, 0.1):
        for y in np.arange(y0-1, yn+2, 0.1):
            if (check_left and is_in((x, y), lst_pixel_ro=lst_pixel_ro)) or \
               ((not check_left) and is_in((x, y), lst_pixel_ro=lst_pixel_ro[::-1])):
                #in_pixels.add((x, y))
                add_nbd_to_set((x, y), set_pixel=in_pixels)

    lst_voxel = [(int(x-0.5), voxel_y, int(y-0.5)) for (x, y) in in_pixels]
    return lst_voxel

def lawson_csv_to_voxel(filename):
    if not os.path.isfile(filename):
        print('File not found: %s' % filename)
        return None
    lst_voxel = []
    lst_voxel_lr = {'r':[], 'l':[]}

    lr = 'r'
    with open(filename) as fin:
        line = fin.readline() # Right
        if line[:5] != 'Right':
            print('1st line is not Right')
        line = fin.readline()
        while line[:3] != 'End':
            lst = [float(value) for value in line.strip().split(',')]
            lst_voxel += pixel_to_voxel(lst, lr=lr)
            lst_voxel_lr[lr] += pixel_to_voxel(lst, lr=lr)
            line = fin.readline()
            if line[:4] == 'Left':
                line = fin.readline()
                lr = 'l'
    return lst_voxel, lst_voxel_lr

lst_voxel, lst_voxel_lr = lawson_csv_to_voxel(fn_csv)
img = nib.load(fn_ref)
if img.affine[0][0] > 0:
    print 'WARNING: Check image orientation'
hdr = img.header
dat_out = np.zeros(img.shape, dtype=np.int8)
hdr.set_data_dtype(np.int8)
for voxel in lst_voxel:
    dat_out[voxel] = 1

if filter_low_intensity:
    dat_ref = img.get_data()
    ary_intensity = dat_ref[dat_out > 0]
    ary_intensity.sort()
    q1 = ary_intensity[ary_intensity.shape[0]/4]
    q3 = ary_intensity[ary_intensity.shape[0]*3/4]
    iqr = q3 - q1
    of1 = q1 - 1.5*iqr
    dat_out[dat_ref < of1] = 0

img_out = nib.Nifti1Image(dat_out, img.affine, hdr)
nib.save(img_out, fn_out)

for lr in ['r', 'l']:
    dat_out_lr = np.zeros(img.shape, dtype=np.int8)
    for voxel in lst_voxel_lr[lr]:
        dat_out_lr[voxel] = 1
    dat_out_lr[dat_out == 0] = 0
    img_out = nib.Nifti1Image(dat_out_lr, img.affine, hdr)
    nib.save(img_out, '%s_%s.nii.gz' % (bn_out, lr))

