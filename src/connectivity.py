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


############ collect and summary the output to dataframe ##############
## version 1
# import pickle
# from tqdm import tqdm
# import numpy as np
# import pandas as pd
# from scipy.spatial.distance import cdist
# import os
#
# if not os.path.exists(
#         r'D:\2022\whole_brain_con\final_data\connectivity_results\Allen_Janelia_re10um_single_cell_level_result_z_unfilpped.pickle'):
#     src = r'D:\2022\whole_brain_con\final_data\connectivity_results\Allen_Janelia_re10um_single_cell_level_result_z_unfilpped'
#     file_list = os.listdir(src)
#     not_list = []
#     result_dict = {}
#     for file in tqdm(file_list):
#         path_file = src + '\\' + file
#         with open(path_file, 'rb') as f:
#             tmp = pickle.load(f)
#         if not len(tmp) == 2838:
#             not_list.append(path_file)
#         cur_result = []
#         for tmp_i in tmp:
#             cur_result.append(
#                 [tmp_i[0].mean(),
#                  np.median(tmp_i[0]),
#                  tmp_i[0][tmp_i[0] < np.percentile(tmp_i[0], 25)].mean(),
#                  tmp_i[0][tmp_i[0] < 50].mean() if len(tmp_i[0][tmp_i[0] < 50]) else 0,
#                  tmp_i[0][tmp_i[0] < 40].mean() if len(tmp_i[0][tmp_i[0] < 40]) else 0,
#                  tmp_i[0][tmp_i[0] < 30].mean() if len(tmp_i[0][tmp_i[0] < 30]) else 0,
#                  tmp_i[0][tmp_i[0] < 20].mean() if len(tmp_i[0][tmp_i[0] < 20]) else 0,
#                  tmp_i[0][tmp_i[0] < 10].mean() if len(tmp_i[0][tmp_i[0] < 10]) else 0,
#                  tmp_i[1],
#                  tmp_i[2]])
#         result_dict[file] = cur_result
#     with open(
#             r'D:\2022\whole_brain_con\final_data\connectivity_results\Allen_Janelia_re10um_single_cell_level_result_z_unfilpped.pickle',
#             'wb') as f:
#         pickle.dump(result_dict, f)
# else:
#     with open(
#             r'D:\2022\whole_brain_con\final_data\connectivity_results\Allen_Janelia_re10um_single_cell_level_result_z_unfilpped.pickle',
#             'rb') as f:
#         result_dict = pickle.load(f)
#
# sele_index = [i.split('.')[0] for i in result_dict.keys()]
# tmp_list1 = result_dict[sele_index[0] + '.pickle']
# sele_column = [j.split('.')[0] for j in list(set([i[-1] for i in tmp_list1] + [tmp_list1[0][-2]]))]
# type_dict = {0: 'mean', 1: 'median', 2: 'mean-per25', 3: 'th50', 4: 'th40', 5:'th30', 6:'th20', 7:'th10'}
# result_dataframe_dict = {}
# for k, v in type_dict.items():
#     tmp_result_df = None
#     for axon_id, tmp_dataset in result_dict.items():
#         values_list = []
#         column_list = []
#         for tmp_data in tmp_dataset:
#             values_list.append(tmp_data[k])
#             column_list.append(tmp_data[-1].split('.')[0])
#         if axon_id.split('.')[0] not in column_list:
#             values_list.append(np.average(values_list))
#             column_list.append(axon_id.split('.')[0])
#         tmp_df = pd.DataFrame(np.array(values_list).reshape(1, -1), columns=column_list, index=[axon_id])[sele_column]
#         if tmp_result_df is None:
#             tmp_result_df = tmp_df
#         else:
#             tmp_result_df = pd.concat([tmp_result_df, tmp_df], axis=0)
#     result_dataframe_dict[v] = tmp_result_df
#
# with open(
#         r'D:\2022\whole_brain_con\final_data\connectivity_results\Allen_Janelia_re10um_single_cell_level_summaryz_unfilpped.pickle',
#         'wb') as f:
#     pickle.dump(result_dataframe_dict, f)
