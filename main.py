# Start date = April 14th
import itertools

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

mobidb_features_lst = []
mobidb_predictors_cont_fra_dict = {}
cont_fra_temp_lst = []

disease_cont_fra_temp_lst = []
disease_predictors_cont_fra_dict = {}


def drawplot(plot_input_lst, bins, is_dense, x_label, y_label, png_file_name, subdirectory):
    plt.hist(plot_input_lst, bins=bins,
             density=is_dense)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    sub_directory = subdirectory  # homosapiens or disease
    file_format = '.png'
    plt.savefig('plots/' + sub_directory + '/' + png_file_name + file_format)
    plt.show()


def compare_plot(first_lst, second_lst, bins, is_dense, first_label, second_label, x_label, y_label, png_file_name):
    plt.hist([first_lst, second_lst], bins=bins, density=is_dense, label=[first_label, second_label])
    plt.legend(loc='upper right')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    file_format = '.png'
    plt.savefig('plots/compare/' + png_file_name + '_compare' + file_format)
    plt.show()


def matrix_maker_nan(input_df, num_3rd_dim):
    matrix_2d = (input_df.to_numpy() <= 1.) * input_df.to_numpy()  # <=1 cuz it indicates disorder percentage(c_f)
    matrix_3d = np.full((matrix_2d.shape[0], matrix_2d.shape[1], num_3rd_dim + 1), np.nan)
    for i in range(matrix_2d.shape[0]):
        for j in range(matrix_2d.shape[1]):
            if matrix_2d[i, j] != 0:
                k = int(round(matrix_2d[i, j] * num_3rd_dim))
                matrix_3d[i, j, k] = 1
    # Replace NaN with zeros for rows containing at least one value
    for i in range(matrix_3d.shape[0]):
        for j in range(matrix_3d.shape[1]):
            if 1.0 in matrix_3d[i][j]:
                matrix_3d[i][j][np.isnan(matrix_3d[i][j])] = 0
    # sum of Pr.s with same content_fraction for each feature
    matrix_3d_sum = np.nansum(matrix_3d, axis=0)
    matrix_3d_sum = matrix_3d_sum / matrix_3d_sum.max(axis=1)[:, None]
    return matrix_2d, matrix_3d, matrix_3d_sum


def matrix_maker_zeros(input_df, num_3rd_dim):
    matrix_2d = (input_df.to_numpy() <= 1.) * input_df.to_numpy()
    matrix_3d = np.zeros(matrix_2d.shape[0], matrix_2d[1], num_3rd_dim + 1)
    for i in range(matrix_2d.shape[0]):
        for j in range(matrix_2d.shape[1]):
            if matrix_2d[i, j] != 0:
                k = int(round(matrix_2d[i, j] * num_3rd_dim))
                matrix_3d[i, j, k] = 1
    matrix_3d_sum = np.sum(matrix_3d, axis=0)
    matrix_3d_sum = matrix_3d_sum / matrix_3d_sum.max(axis=1)[:, None]
    return matrix_2d, matrix_3d, matrix_3d_sum


def sum_df_generator(input_sum_matrix, index_list):
    sum_df = pd.DataFrame(input_sum_matrix,
                          columns=['0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'],
                          index=index_list)
    return sum_df


def draw_heatmaps(sum_df_1, title_1, sum_df_2, title_2, ylabel, xlabel, saving_rout):
    sns.set()
    fig, ax = plt.subplots(2, 1, figsize=(24, 12))
    for i, d in enumerate([sum_df_1, sum_df_2]):
        labels = d.applymap(lambda v: str(v) if v == d.values.max() else '')
        sns.heatmap(d,
                    cmap="viridis",  # sequential colormap
                    annot_kws={'fontsize': 11},
                    fmt='',
                    square=True,
                    # vmax=1,
                    # vmin=0,
                    linewidth=0.01,
                    linecolor="#222",
                    ax=ax[i],
                    )
    ax[0].set_title(title_1)
    ax[1].set_title(title_2)
    ax[0].set_ylabel(ylabel)
    ax[1].set_ylabel(ylabel)
    ax[0].set_xlabel(xlabel)
    ax[1].set_xlabel(xlabel)
    plt.tight_layout()
    plt.savefig(saving_rout, dpi=120)
    plt.show()


# if __name__ == '__main__':
## Files import and manipulation
mobidb_original_df = pd.read_csv('data/mobidb_result.tsv', sep='\t')
mobidb_transposed_df = mobidb_original_df.pivot_table(index=['acc'], columns=['feature'],
                                                      values='content_fraction').fillna(0)
mobidb_transposed_df = mobidb_transposed_df.reset_index()  # added idx nums manually,ACCs recognized separate column
mobidb_transposed_df.to_csv(r'data/mobidb_transposed_df.csv', index=True)
mobidb_features_lst = mobidb_transposed_df.columns.str.split(',').tolist()  # this also contains the 'acc' column
mobidb_features_lst = list(itertools.chain(*mobidb_features_lst))  # flat list

disease_acc_df = pd.read_csv('data/diseases.tab', sep='\t')
disease_acc_lst = disease_acc_df['Entry'].to_list()
disease_mobidb_df = mobidb_transposed_df[mobidb_transposed_df['acc'].isin(disease_acc_lst)]

## Matrix
mobidb_matrix, mobidb_3d_matrix, mobidb_3d_matrix_sum = matrix_maker_nan(mobidb_transposed_df.iloc[:, 1:], 10)
disease_matrix, disease_3d_matrix, disease_3d_matrix_sum = matrix_maker_nan(disease_mobidb_df.iloc[:, 1:], 10)

## Sum df for heat map
mobidb_cont_fract_sum_df = sum_df_generator(mobidb_3d_matrix_sum, mobidb_features_lst[1:])
disease_cont_fract_sum_df = sum_df_generator(disease_3d_matrix_sum, mobidb_features_lst[1:])
#TODO: find the difference of two matrices

## heatmaps
draw_heatmaps(sum_df_1=mobidb_cont_fract_sum_df.T, title_1='Homo sapiens', sum_df_2=disease_cont_fract_sum_df.T,
              title_2='NDDs', ylabel='Disorder %', xlabel='mobiDB features', saving_rout='plots/heatmaps/Hmaps.png')

## Dictionary Homo sapiens
for each_feature in mobidb_features_lst:
    cont_fra_temp_lst = mobidb_transposed_df[each_feature].tolist()
    mobidb_predictors_cont_fra_dict[each_feature] = cont_fra_temp_lst

## Plot for homosapiens
for each_feature in mobidb_features_lst[
                    1:]:
    drawplot(mobidb_predictors_cont_fra_dict[each_feature], 30, True, each_feature + '_homosapiens',
             "Protein Count(relative)",
             each_feature, 'homosapiens')

## Dictionary Disease
for each_feature in mobidb_features_lst:
    disease_cont_fra_temp_lst = disease_mobidb_df[each_feature].tolist()
    disease_predictors_cont_fra_dict[each_feature] = disease_cont_fra_temp_lst

## plot for diseases
for each_feature in mobidb_features_lst[1:]:
    drawplot(disease_predictors_cont_fra_dict[each_feature], 30, False, each_feature + '_Disease', 'Protein count',
             each_feature + '_disease', 'disease')

## comparative histogram (homosapiens Vs. disease)
for each_feature in mobidb_features_lst[
                    1:]:
    compare_plot(mobidb_predictors_cont_fra_dict[each_feature], disease_predictors_cont_fra_dict[each_feature], 30,
                 True, x_label=each_feature + '_comparison', y_label='proteins count(relative)',
                 png_file_name=each_feature,
                 first_label='all_disordered_Pr.s', second_label='disease_Pr.s')
