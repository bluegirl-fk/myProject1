# Start date = April 14th
import itertools

import numpy as np
# import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

mobidb_features_lst = []
mobidb_predictors_cont_fra_dict = {}
cont_fra_temp_lst = []

ndd_cont_fra_temp_lst = []
ndd_predictors_cont_fra_dict = {}


def drawplot(plot_input_lst, bins, is_dense, x_label, y_label, png_file_name, subdirectory):
    plt.hist(plot_input_lst, bins=bins,
             density=is_dense)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    sub_directory = subdirectory  # homosapiens or ndd
    file_format = '.png'
    plt.savefig('plots/' + sub_directory + '/' + png_file_name + file_format)
    plt.show()
    return


def compare_plot(first_lst, second_lst, bins, is_dense, first_label, second_label, x_label, y_label, png_file_name):
    plt.hist([first_lst, second_lst], bins=bins, density=is_dense, label=[first_label, second_label])
    plt.legend(loc='upper right')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    file_format = '.png'
    plt.savefig('plots/' + png_file_name + '_compare' + file_format)
    plt.show()
    return

def matrix_maker_nan(input_df, num_3rd_dim):
    matrix_2d = (input_df.to_numpy() <= 1.) * input_df.to_numpy()  # <=1 cuz it indicates disorder percentage(c_f)
    matrix_3d = np.full((matrix_2d.shape[0], matrix_2d.shape[1], num_3rd_dim + 1), np.nan)
    for i in range(matrix_2d.shape[0]):
        for j in range(matrix_2d.shape[1]):
            if matrix_2d[i, j] != 0:
                k = int(round(matrix_2d[i, j] * num_3rd_dim))
                matrix_3d[i, j, k] = 1

    #
    # Replace NaN with zeros for rows containing at least one value
    # for i in range(matrix_3d.shape[0]):
    #     for j in range(matrix_3d.shape[1]):
    #         if 1.0 in matrix_3d[i][j]:
    #             matrix_3d[i][j][np.isnan(matrix_3d[i][j])] = 0

    # sum of Pr.s with same content_fraction for each feature
    matrix_3d_sum = np.nansum(matrix_3d, axis=0)
    matrix_3d_sum_normalized = matrix_3d_sum / matrix_3d_sum.max(axis=1)[:, None]
    return matrix_2d, matrix_3d, matrix_3d_sum, matrix_3d_sum_normalized


def matrix_maker_zeros(input_df, num_3rd_dim):
    matrix_2d = (input_df.to_numpy() <= 1.) * input_df.to_numpy()
    matrix_3d = np.zeros((matrix_2d.shape[0], matrix_2d.shape[1], num_3rd_dim + 1))
    for i in range(matrix_2d.shape[0]):
        for j in range(matrix_2d.shape[1]):
            # if matrix_2d[i, j] != 0:
            k = int(round(matrix_2d[i, j] * num_3rd_dim))
            matrix_3d[i, j, k] = 1
    matrix_3d_sum = np.sum(matrix_3d, axis=0)
    matrix_3d_sum_normalized = matrix_3d_sum / matrix_3d_sum.max(axis=1)[:, None]
    return matrix_2d, matrix_3d, matrix_3d_sum, matrix_3d_sum_normalized


def sum_df_generator(input_sum_matrix):
    sum_df = pd.DataFrame(input_sum_matrix,
                          columns=['0', ' ', '20', ' ', '40', ' ', '60', ' ', '80', ' ', '100'],
                          index=mobidb_features_lst[1:])
    return sum_df


# def draw_heatmaps(data, titles, saving_rout):  # from stackabuse.com/ultimate-guide-to-heatmaps-in-seaborn-with-python/
#     sns.set()
#     fig, axes = plt.subplots(len(data), 1, figsize=(12 * len(data), 12))
#     for i, (ax, d, t) in enumerate(zip(axes.reshape(-1), data, titles)):
#         sb = sns.heatmap(d,
#                          cmap="viridis",  # sequential colormap
#                          annot_kws={'fontsize': 11},
#                          fmt='',
#                          square=True,
#                          # vmax=1,
#                          # vmin=0,
#                          linewidth=0.01,
#                          linecolor="#222",
#                          ax=ax,
#                          vmin=-1.0, vmax=1.0
#                          )
#         ax.set_title(t)
#         # ax.set_ylabel(ylabel)
#         if i < (len(data) - 1):
#             sb.set(xticklabels=[])
#             sb.set(xlabel=None)
#     plt.tight_layout()
#     plt.savefig(saving_rout, dpi=120)
#     plt.show()

    return


# if __name__ == '__main__':
## Files import and manipulation
mobidb_original_df = pd.read_csv('data/mobidb_result.tsv', sep='\t')
mobidb_transposed_df = mobidb_original_df.pivot_table(index=['acc'], columns=['feature'],
                                                      values='content_fraction').fillna(0)
mobidb_transposed_df = mobidb_transposed_df.reset_index()  # added idx nums manually,ACCs recognized separate column
mobidb_transposed_df.to_csv(r'data/mobidb_transposed_df.csv', index=True)
mobidb_features_lst = mobidb_transposed_df.columns.str.split(',').tolist()  # this also contains the 'acc' column
mobidb_features_lst = list(itertools.chain(*mobidb_features_lst))  # flat list

ndd_acc_df = pd.read_csv('data/allUniqueEntry.tab', sep='\t')
ndd_acc_lst = ndd_acc_df['Entry'].to_list()
ndd_mobidb_df = mobidb_transposed_df[mobidb_transposed_df['acc'].isin(ndd_acc_lst)]

## Matrix
# with nan
_, mobidb_3d_matrix_nan, mobidb_3d_matrix_nan_sum, mobidb_3d_matrix_nan_sum_norm = matrix_maker_nan(
    mobidb_transposed_df.iloc[:, 1:], 10)
_, ndd_3d_matrix_nan, ndd_3d_matrix_nan_sum, ndd_3d_matrix_nan_sum_norm = matrix_maker_nan(ndd_mobidb_df.iloc[:, 1:],
                                                                                           10)

## Sum_norm df for heat map
mobidb_cont_fract_sum_norm_df = sum_df_generator(mobidb_3d_matrix_nan_sum_norm)
ndd_cont_fract_sum_norm_df = sum_df_generator(ndd_3d_matrix_nan_sum_norm)

## Difference of the sum arrays(with nan)
sum_difference_matrix_nan_norm = mobidb_3d_matrix_nan_sum_norm - ndd_3d_matrix_nan_sum_norm
sum_difference_df_nan_norm = sum_df_generator(sum_difference_matrix_nan_norm)

## heatmaps
# draw_heatmaps([mobidb_cont_fract_sum_norm_df.T, ndd_cont_fract_sum_norm_df.T, sum_difference_df_nan_norm.T],
#               ['Homo sapiens', 'Attention Deficit Hyperactivity Disorder(ADHD)', 'Difference (Homo sapiens - ADHD)'],
#               saving_rout='plots/heatmaps/Hmaps-ADHD.png')

## Sum (not normalized) df for stacked histogram
mobidb_3d_matrix_nan_sum_df = sum_df_generator(mobidb_3d_matrix_nan_sum)
ndd_3d_matrix_nan_sum_df = sum_df_generator(ndd_3d_matrix_nan_sum)

# # data for stacked histogram of mobiDB and NDD separately x1_0_perc_lst = mobidb_3d_matrix_nan_sum_df['0'].to_list(
# ) x1_20_perc_lst = mobidb_3d_matrix_nan_sum_df['20'].to_list() x1_40_perc_lst = mobidb_3d_matrix_nan_sum_df[
# '40'].to_list() x1_100_perc_lst = mobidb_3d_matrix_nan_sum_df['100'].to_list() plt.figure() plt.hist([
# x1_0_perc_lst, x1_20_perc_lst, x1_40_perc_lst, x1_100_perc_lst],mobidb_features_lst[1:], bins=20, stacked=True,
# density=True) plt.show() apply log() separated 3d_matrix_sum and sum_normalized variables

# import sys
#
# sys.exit(0)

## Dictionary Homo sapiens
for each_feature in mobidb_features_lst:
    cont_fra_temp_lst = mobidb_transposed_df[each_feature].tolist()
    mobidb_predictors_cont_fra_dict[each_feature] = cont_fra_temp_lst

## Plot for homosapiens
# for each_feature in mobidb_features_lst[
#                     1:]:
#     drawplot(mobidb_predictors_cont_fra_dict[each_feature], 30, True, each_feature + '_homosapiens',
#              "Protein Count(relative)",
#              each_feature, 'homosapiens')

## Dictionary Disease
for each_feature in mobidb_features_lst:
    ndd_cont_fra_temp_lst = ndd_mobidb_df[each_feature].tolist()
    ndd_predictors_cont_fra_dict[each_feature] = ndd_cont_fra_temp_lst

## plot for ndds
for each_feature in mobidb_features_lst[1:]:
    drawplot(ndd_predictors_cont_fra_dict[each_feature], 30, False, each_feature + '_NDD', 'Protein count',
             each_feature + '_NDD', 'NDD-all-hist/NDD')

## comparative histogram (homosapiens Vs. ndd)
for each_feature in mobidb_features_lst[
                    1:]:
    compare_plot(mobidb_predictors_cont_fra_dict[each_feature], ndd_predictors_cont_fra_dict[each_feature], 30,
                 True, x_label=each_feature + '_comparison', y_label='proteins count(relative)',
                 png_file_name='/NDD-all-hist/compare/'+each_feature,
                 first_label='Homo sapiens Pr.s', second_label='NDD Pr.s')
