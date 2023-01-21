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


def cal_connectivity(swc1, swc2, swc2_soma, swc_id1, swc_id2):
    try:
        # print(swc_id1, swc_id2)
        axon1 = swc1.loc[swc1.type.isin([2]), ['x', 'y', 'z']].values
        den2 = swc2.loc[swc2.type.isin([3, 4]), ['x', 'y', 'z']].values
        soma2 = swc2_soma.iloc[0].values
        if len(axon1) > 1 and len(den2) > 1:
            distance_matrix = cdist(axon1, den2)  # scipy
            distance = distance_matrix.min(axis=1)

            distance_1 = (axon1 - soma2)
            distance_2 = distance_1*distance_1
            distance_n2s = np.sqrt(distance_2.sum(axis=1))

            return distance, distance_n2s, swc_id1, swc_id2
        else:
            return -1, swc_id1, swc_id2

    except:
        print('fail at: ', swc_id1, ' ', swc_id2)
        return 0, swc_id1, swc_id2


# def cal_node2soma(swc1, swc2, swc_id1, swc_id2):
#     try:
#         # print(swc_id1, swc_id2)
#         axon1 = swc1.loc[swc1.type.isin([2]), ['x', 'y', 'z']].values
#         soma2 = swc2.iloc[0].values
#         if len(axon1) > 1:
#             distance_1 = (axon1 - soma2)
#             distance_2 = distance_1*distance_1
#             distance_n2s = np.sqrt(distance_2.sum(axis=1))
#             return distance_n2s, swc_id1, swc_id2
#         else:
#             print(swc_id1, ' ', swc_id2, ' error')
#             return -1, swc_id1, swc_id2
#
#     except:
#         print('fail at: ', swc_id1, ' ', swc_id2)
#         return 0, swc_id1, swc_id2


if __name__ == "__main__":
    src1 = r'D:\datasets\Allen_Janelia\Allen_Jan_re10um'
    swc_list1 = os.listdir(src1)
    swc_path_list1 = [src1 + '\\' + i for i in swc_list1]
    swc_dict1 = {}
    for cur_swc_path in tqdm(swc_path_list1):
        swc_dict1[cur_swc_path.split('\\')[-1]] = readSWC(cur_swc_path)

    soma_df = pd.read_csv('./final_data/3002_cell_soma.csv', sep=',', index_col=0)
    denRadius = pd.read_csv('./final_data/denRadius.csv', sep=',', index_col=0)
    # denRadius

    cores = int(multiprocessing.cpu_count() * 0.8)  # multiprocessing.cpu_count()
    print("cores: ", cores)

    # param_list = []
    for swc_i, dfi in swc_dict1.items():
        if os.path.exists(r'D:\datasets\result' + '\\' + swc_i.split('.')[0] + '.pickle'):
            continue
        print('running at: ', swc_i)
        pool = multiprocessing.Pool(processes=cores)
        res = []
        for swc_j, dfj in (swc_dict1.items()):
            if swc_i == swc_j:
                continue
            # param_list.append([df1, df2, swc_i, swc_j])
            swcj_soma = soma_df.loc[[swc_j], ['x', 'y', 'z']]
            res.append(pool.apply_async(cal_connectivity, (dfi, dfj, swcj_soma, swc_i, swc_j)))
            # _, _, _ = cal_connectivity(df1, df2, swc_i, swc_j)

        # for param in param_list:
        #     res.append(pool.apply_async(cal_connectivity, (param[0], param[1], param[2], param[3])))

        # result_df = pd.DataFrame(np.zeros((len(swc_list), len(swc_list))), index=list(swc_dict.keys()),
        #                          columns=list(swc_dict.keys()))

        cur_result = []
        for tmp_res in res:
            tmp_i = tmp_res.get()
            if type(tmp_i[0]) == np.ndarray:
                axon_id = tmp_i[2]
                den_id = tmp_i[3]
                a2d_dis = tmp_i[0]
                a2s_dis = tmp_i[1]
                cur_radius = denRadius.loc[den_id.split('.')[0]].values[0]
                cur_index = np.argwhere(a2s_dis < cur_radius)

                cur_result.append(
                    [round(a2d_dis.mean(), 3),
                     round(a2d_dis[cur_index].mean(), 3) if len(cur_index) > 0 else -1,
                     len(cur_index),
                     axon_id,
                     den_id])

        with open(r'D:\datasets\result' + '\\' + swc_i.split('.')[0] + '.pickle',
                  'wb') as f:
            pickle.dump(cur_result, f)

        print('finish : ', swc_i)
        #
        # # result_df.to_csv(src + '.csv', sep=',')
        #
        pool.close()  # 关闭进程池，表示不能在往进程池中添加进程
        pool.join()  # 等待进程池中的所有进程执行完毕，必须在close()之后调用


