# Start date = April 14th
import itertools

import matplotlib.pyplot
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

mobidb_features_lst = []
mobidb_predictors_cont_fra_dict = {}
cont_fra_temp_lst = []

ndd_cont_fra_temp_lst = []
ndd_predictors_cont_fra_dict = {}


def drawplot(plot_input_lst, yscale, bins, is_dense, x_label, y_label, png_file_name):
    plt.hist(plot_input_lst, bins=bins,
             density=is_dense)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    file_format = '.png'
    plt.yscale(yscale)
    plt.savefig(png_file_name + file_format)
    plt.show()
    return


def draw_barplot(figsize_a, figsize_b, xlabel, ylabel, data, xticklabel, yscale, save_rout):  # input is DF, not list
    plt.figure(figsize=(int(figsize_a), int(figsize_b)))  # bigger figsize to have xticklabels shown
    sns.set_style("ticks")
    g = sns.barplot(x=xlabel, y=ylabel, data=data)
    # sns.despine(trim=True, offset=2)
    g.set_xticklabels(xticklabel)
    sns.color_palette("pastel")
    plt.yscale(yscale)
    plt.tight_layout()
    plt.savefig(save_rout)
    plt.show()
    return


def compare_plot(first_lst, second_lst, yscale, bins, is_dense, first_label, second_label, x_label, y_label,
                 png_file_name):
    plt.hist([first_lst, second_lst], bins=bins, density=is_dense, label=[first_label, second_label])
    plt.legend(loc='upper right')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    file_format = '.png'
    plt.yscale(yscale)
    plt.savefig(png_file_name + '_compare' + file_format)
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


def draw_heatmaps(data, titles, saving_rout):  # www.stackabuse.com/ultimate-guide-to-heatmaps-in-seaborn-with-python/
    sns.set()
    fig, axes = plt.subplots(len(data), 1, figsize=(12 * len(data), 12))
    for i, (ax, d, t) in enumerate(zip(axes.reshape(-1), data, titles)):
        sb = sns.heatmap(d,
                         cmap="viridis",  # sequential colormap
                         annot=False,
                         annot_kws={'fontsize': 6},
                         fmt='',
                         square=True,
                         # vmax=1,
                         # vmin=0,
                         linewidth=0.01,
                         linecolor="#222",
                         ax=ax,
                         vmin=-1.0, vmax=1.0
                         )
        ax.set_title(t)
        # ax.set_ylabel(ylabel)
        if i < (len(data) - 1):
            sb.set(xticklabels=[])
            sb.set(xlabel=None)
    plt.tight_layout()
    plt.savefig(saving_rout, dpi=120)
    plt.show()

    return


# if __name__ == '__main__':
### Files import and modify
mobidb_original_df = pd.read_csv('data/mobidb_result.tsv', sep='\t')
## for content fraction
mobidb_pivot_contf_df = mobidb_original_df.pivot_table(index=['acc'], columns=['feature'],
                                                       values='content_fraction').fillna(0)
mobidb_pivot_contf_df = mobidb_pivot_contf_df.reset_index()  # added idx nums manually,ACCs recognized separate column
mobidb_pivot_contf_df.to_csv(r'data/mobidb_pivot_contf_df.csv', index=True)
mobidb_features_lst = mobidb_pivot_contf_df.columns.str.split(',').tolist()  # this also contains the 'acc' column
mobidb_features_lst = list(itertools.chain(*mobidb_features_lst))  # flat list

ndd_acc_df = pd.read_csv('data/allUniqueEntry.tab', sep='\t')
ndd_acc_lst = ndd_acc_df['Entry'].to_list()
ndd_mobidb_df = mobidb_pivot_contf_df[mobidb_pivot_contf_df['acc'].isin(ndd_acc_lst)]
## for Length
mobidb_pivot_length_df = mobidb_original_df.pivot_table(index=['acc'], columns=['feature'],
                                                        values='length').fillna(0)  # (75052, 78)

## Matrix
# with nan
_, mobidb_3d_matrix_nan, mobidb_3d_matrix_nan_sum, mobidb_3d_matrix_nan_sum_norm = matrix_maker_nan(
    mobidb_pivot_contf_df.iloc[:, 1:], 10)
_, ndd_3d_matrix_nan, ndd_3d_matrix_nan_sum, ndd_3d_matrix_nan_sum_norm = matrix_maker_nan(ndd_mobidb_df.iloc[:, 1:],
                                                                                           10)

# # Add the length statistics here. Use vstack or hstack
# # mobidb_3d_matrix_nan_sum
# # mobidb_3d_matrix_nan_sum_norm
# # ndd_3d_matrix_nan_sum
# # ndd_3d_matrix_nan_sum_norm
#
# #ax.hist(dataset_len, bins=np.arange(0, 1000, 10))
#
# # for length
# matrix_2d_length = (mobidb_pivot_length_df.to_numpy() <= 1.) * mobidb_pivot_length_df.to_numpy()
# matrix_3d_length = np.zeros((matrix_2d_length.shape[0], matrix_2d_length.shape[1], 10 + 1))
# for i in range(matrix_2d_length.shape[0]):
#     for j in range(matrix_2d_length.shape[1]):
#         # if matrix_2d[i, j] != 0:
#         k = int(round(matrix_2d_length[i, j] * 10))
#         matrix_3d_length[i, j, k] = 1
# matrix_3d_sum = np.sum(matrix_3d_length, axis=0)
# matrix_3d_sum_normalized = matrix_3d_sum / matrix_3d_sum.max(axis=1)[:, None]


## columns sum of matrix_3d_sum df to get prot count per feature (for histogram based on distribution of heatmap)
mobidb_columns_sum_df = pd.DataFrame([mobidb_3d_matrix_nan_sum.T.sum(axis=0)], columns=mobidb_features_lst[1:],
                                     index=['Proteins count'])
ndd_columns_sum_df = pd.DataFrame([ndd_3d_matrix_nan_sum.T.sum(axis=0)], columns=mobidb_features_lst[1:],
                                  index=['Proteins count'])
mobidb_columns_sum_df = mobidb_columns_sum_df.T.reset_index()
ndd_columns_sum_df = ndd_columns_sum_df.T.reset_index()
mobidb_columns_sum_df.columns = ['Features', 'Protein count']
ndd_columns_sum_df.columns = ['Features', 'Protein count']
mobidb_cols_sum_lst = mobidb_columns_sum_df['Protein count']
mobidb_cols_sum_lst = [int(x) for x in mobidb_cols_sum_lst]
ndd_cols_sum_lst = ndd_columns_sum_df['Protein count']
ndd_cols_sum_lst = [int(x) for x in ndd_cols_sum_lst]
# ## Gene4denovo  (delete acc duplicates)
# gene4dn_all_annotations_df = pd.read_csv('data/gene4denovo/All_De_novo_mutations_and_annotations_1.2.txt',
#                                          sep='\t', encoding='cp1252', low_memory=False)  # (670082, 155)
# # filter this columns : exonic,ENSG00000115020,-,nonsynonymous SNV,
# # ENSG00000115020:ENST00000452564:exon19:c.2929T>A:p.S977T,
# # ENSG00000115020:ENST00000264380:exon20:c.3097T>A:p.S1033T,-,-,-,-,-,-,-,-,-,-,-,-,-,
# genes4dn_orig_df = pd.read_csv('data/gene4denovo/genes4dn.txt', sep='\t')  # (8271, 13)
# genes4dn_acc_df = pd.read_csv('data/uniprot-gene4dn-acc.tab', sep='\t')  # (8039, 7)
# genes4dn_acc_merge_df = pd.merge(genes4dn_orig_df, genes4dn_acc_df, on='geneslist')  # (48060, 19)

## sum dataframes
mobidb_cont_fract_sum_norm_df = sum_df_generator(mobidb_3d_matrix_nan_sum_norm)
ndd_cont_fract_sum_norm_df = sum_df_generator(ndd_3d_matrix_nan_sum_norm)

## Difference of the sum arrays(with nan)
sum_difference_matrix_nan_norm = mobidb_3d_matrix_nan_sum_norm - ndd_3d_matrix_nan_sum_norm
sum_difference_df_nan_norm = sum_df_generator(sum_difference_matrix_nan_norm)

## heatmaps
# draw_heatmaps([mobidb_cont_fract_sum_norm_df.T, ndd_cont_fract_sum_norm_df.T, sum_difference_df_nan_norm.T],
#               ['Homo sapiens', 'NDDs', 'Difference (Homo sapiens - NDDs)'],
#               saving_rout='plots/heatmaps/Hmaps1.png')
#
# mobidb_cont_fract_sum_norm_df.index = mobidb_cont_fract_sum_norm_df.index.set_names(['Features'])
# ndd_cont_fract_sum_norm_df.index = ndd_cont_fract_sum_norm_df.index.set_names(['Features'])
# merged_mobidb_hmap_df = pd.merge(mobidb_columns_sum_df, mobidb_cont_fract_sum_norm_df, on='Features').set_index(
#     'Features')
# merged_ndd_hmap_df = pd.merge(ndd_columns_sum_df, ndd_cont_fract_sum_norm_df, on='Features').set_index('Features')
#
# draw_heatmaps([merged_mobidb_hmap_df.T, merged_ndd_hmap_df.T, sum_difference_df_nan_norm.T],
#               ['Homo sapiens', 'NDDs', 'Difference (Homo sapiens - NDDs)'],
#               saving_rout='plots/heatmaps/Hmap_with_sum.png')

## distribution heatmap plot
draw_barplot(figsize_a='40', figsize_b='20', xlabel='Features', ylabel='Protein count', data=mobidb_columns_sum_df,
             xticklabel=mobidb_cols_sum_lst, yscale='log', save_rout='plots/log/hist-hmaps-distribution/mobidb-log.png')
draw_barplot(figsize_a='40', figsize_b='20', xlabel='Features', ylabel='Protein count', data=ndd_columns_sum_df,
             xticklabel=ndd_cols_sum_lst, yscale='log', save_rout='plots/log/hist-hmaps-distribution/ndd-log.png')

# Protein count
import sys

sys.exit(0)
## Dictionary Homo sapiens
for each_feature in mobidb_features_lst:
    cont_fra_temp_lst = mobidb_pivot_contf_df[each_feature].tolist()
    mobidb_predictors_cont_fra_dict[each_feature] = cont_fra_temp_lst

## Plot for homosapiens
for each_feature in mobidb_features_lst[
                    1:]:
    drawplot(mobidb_predictors_cont_fra_dict[each_feature], 'log', 30, True, each_feature + '_homosapiens',
             "Protein Count(relative)",
             'plots/log/hist-all-homo sapiens/' + each_feature)

## Dictionary Disease
for each_feature in mobidb_features_lst:
    ndd_cont_fra_temp_lst = ndd_mobidb_df[each_feature].tolist()
    ndd_predictors_cont_fra_dict[each_feature] = ndd_cont_fra_temp_lst

## plot for ndds
for each_feature in mobidb_features_lst[1:]:
    drawplot(ndd_predictors_cont_fra_dict[each_feature], 'log', 30, False, each_feature + '_NDD', 'Protein count',
             'plots/log/hist-all-NDD/SCZ/' + each_feature + '_SCZ')

# comparative histogram (homosapiens Vs. ndd)
for each_feature in mobidb_features_lst[1:]:
    compare_plot(first_lst=mobidb_predictors_cont_fra_dict[each_feature],
                 second_lst=ndd_predictors_cont_fra_dict[each_feature], yscale='log', bins=30, is_dense=False,
                 x_label=each_feature + '_comparison', y_label='proteins count(relative)',
                 first_label='Homo sapiens Pr.s', second_label='NDD Pr.s',
                 png_file_name='plots/log/hist-comparison' '-homoS-NDD/' + each_feature)
