import os
import pandas as pd
import numpy as np


def readSWC(swc_path, mode='simple', scale=1.0):
    n_skip = 0
    with open(swc_path, "r") as f:
        for line in f.readlines():
            line = line.strip()
            if line.startswith("#"):
                n_skip += 1
            else:
                break

    names = ["##n", "type", "x", "y", "z", "r", "parent"]
    used_cols = [0, 1, 2, 3, 4, 5, 6]

    if mode == 'simple':
        pass

    elif mode == 'arbor':
        names = ["##n", "type", "x", "y", "z", "r", "parent", 'arbor_id']
        used_cols = [0, 1, 2, 3, 4, 5, 6, 7]

    elif mode == 'bouton':
        names = ["##n", "type", "x", "y", "z", "r", "parent", 'x_ccf', 'y_ccf', 'z_ccf', 'arbor_id']
        used_cols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    df = pd.read_csv(swc_path, index_col=0, skiprows=n_skip, sep=" ",
                     usecols=used_cols,
                     names=names
                     )
    df['x'] = df['x'] * scale
    df['y'] = df['y'] * scale
    df['z'] = df['z'] * scale

    return df


src1 = r'D:\datasets\3um_SEU_dataset_z_flipped'
src2 = r'D:\datasets\1741_original_s25'
dest = r'D:\datasets\3um_SEU_dataset_z_unflipped'

swc_list = os.listdir(src1)

flip_list = []
for swc in swc_list:
    path1 = src1 + '\\' + swc
    path2 = src2 + '\\' + swc
    if not os.path.exists(path2):
        continue
    df1 = readSWC(path1, scale=0.04)
    df2 = readSWC(path2)
    z1 = df1.loc[(df1.parent == -1), 'z'].values
    z2 = df2.loc[(df2.parent == -1), 'z'].values
    if len(z1) > 1:
        z1 = df1.loc[(df1.type == 1) & (df1.parent == -1), 'z'].values[0]
        z2 = df2.loc[(df2.type == 1) & (df2.parent == -1), 'z'].values[0]
    else:
        z1 = z1[0]
        z2 = z2[0]
    if (z2 > 228) and (z1 < 228):
        flip_list.append(swc)
        df1['z'] = 456 - df1['z']

    df1.to_csv(dest+'\\'+swc, sep=' ')


## batch resample
import os
source = r"D:\datasets\Allen_Janelia\1906_v2"  # swc的文件夹 集合
dest = r"D:\datasets\result"
if not os.path.exists(dest):
    os.mkdir(dest)
bat_path = r"D:\datasets\run_plugin.bat"  # .bat文件地址
#
files = os.listdir(source)
for file in files:
    src = source + "\\" + file
    swc_list = os.listdir(src)
    cur_dest = dest + "\\" + file
    for swc in swc_list:
        if not os.path.exists(cur_dest):
            os.mkdir(cur_dest)
        with open(bat_path, 'w') as OPATH:
            OPATH.writelines(["@echo off",
                              "\n",
                              "D:\\Vaa3D_V3.601_Windows_MSVC_64bit",
                              "\\vaa3d_msvc.exe /x resample_swc /f resample_swc /i ",
                              src + "\\" + swc + " /o ",
                              cur_dest + "\\" + swc + " /p 10"])
        filepath_Aff = bat_path
        os.system(filepath_Aff)


## check multi-soma and duplicates

from tqdm import tqdm
import numpy as np
import pandas as pd
import os

def readSWC(swc_path, mode='simple'):
    n_skip = 0
    with open(swc_path, "r") as f:
        for line in f.readlines():
            line = line.strip()
            if line.startswith("#"):
                n_skip += 1
            else:
                break
    names = ["##n", "type", "x", "y", "z", "r", "parent"]
    used_cols = [0, 1, 2, 3, 4, 5, 6]
    if mode == 'simple':
        pass
    df = pd.read_csv(swc_path, index_col=0, skiprows=n_skip, sep=" ",
                     usecols=used_cols,
                     names=names
                     )

    return df

src=r'D:\datasets\manual_final\ori_1906'

soma_dict = {}
for swc in tqdm(os.listdir(src)):
    df = readSWC(src+'\\'+swc)
    soma = df.loc[(df.parent==-1) & (df.type==1), ['x','y','z']].values
    soma_dict[swc] = soma

dup_list = []
well_list = []
multisoma_list = []
keys_list = list(soma_dict.keys())
for i, k1 in (enumerate(keys_list)):
    for k2 in keys_list[i:]:
        if k1 == k2:
            continue
        if len(soma_dict[k1]) > 1:
            if k1 not in multisoma_list:
                multisoma_list.append(k1)
            continue
        if len(soma_dict[k2]) > 1:
            if k2 not in multisoma_list:
                multisoma_list.append(k2)
            continue
        tmp = soma_dict[k1] - soma_dict[k2]
        dis = np.sqrt((tmp * tmp).sum())
        if dis < 5:
            dup_list.append([k1, k2])


import multiprocessing
import pickle
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import os


def readSWC(swc_path, mode='simple'):
    n_skip = 0
    with open(swc_path, "r") as f:
        for line in f.readlines():
            line = line.strip()
            if line.startswith("#"):
                n_skip += 1
            else:
                break
    names = ["##n", "type", "x", "y", "z", "r", "parent"]
    used_cols = [0, 1, 2, 3, 4, 5, 6]
    if mode == 'simple':
        pass
    df = pd.read_csv(swc_path, index_col=0, skiprows=n_skip, sep=" ",
                     usecols=used_cols,
                     names=names
                     )

    return df
src1=r'D:\datasets\Allen_Janelia\1741_refined'
src2=r'D:\datasets\Allen_Janelia\SEU_V2'
swc_list1 = os.listdir(src1)
swc_list2 = os.listdir(src2)

not_list = []
for swc in swc_list2:
    if swc not in swc_list1:
        not_list.append(swc)



import os
import shutil
src1=r'D:\datasets\Allen_Janelia\SEU_V2\ori'
src2=r'D:\datasets\Allen_Janelia\SEU_V2\rename'
swc_list1=os.listdir(src1)
swc_list2 = os.listdir(src2)
swc_list2 = [i.split('.swc')[0] for i in swc_list2]
new_list = []
for swc in swc_list2:
    swc_split = swc.split('_')
    if swc.count('_stamp'):
        if swc.startswith('pre'):
            new_list.append(swc_split[1]+'_'+swc_split[2]+'.swc')
        else:
            new_list.append(swc_split[0]+'_'+swc_split[1]+'.swc')
    else:
        new_list.append(swc+'.swc')

for i, j in zip(new_list, swc_list1):
    print(i, ' ', j)

for k,i in enumerate(new_list):
    os.rename(src2+'\\'+swc_list1[k], src2+'\\'+i)


# assign each cell to brain-named folder respectively
src=r'D:\datasets\Allen_Janelia\SEU_1892'
swc_list = os.listdir(src)
for swc in swc_list:
    brain_id = swc.split('_')[0]
    brain_folder = src+'\\'+brain_id
    if not os.path.exists(brain_folder):
        os.mkdir(brain_folder)
    shutil.copy(src+'\\'+swc, brain_folder)