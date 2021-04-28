# Start date = April 14th
import itertools

import numpy as np
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
    plt.savefig('graphs/' + sub_directory + '/' + png_file_name + file_format)
    plt.show()


def compare_plot(first_lst, second_lst, bins, is_dense, first_label, second_label, x_label, y_label, png_file_name):
    plt.hist([first_lst, second_lst], bins=bins, density=is_dense, label=[first_label, second_label])
    plt.legend(loc='upper right')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    file_format = '.png'
    plt.savefig('graphs/compare/' + png_file_name + '_compare' + file_format)
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

## 3d matrix for Bayes' theorem
disease_mobidb_matrix = (disease_mobidb_df.iloc[:, 1:].to_numpy() <= 1.) * disease_mobidb_df.iloc[:, 1:].to_numpy()
matrix = np.zeros((disease_mobidb_matrix.shape[0], disease_mobidb_matrix.shape[1], 11))
for i in range(disease_mobidb_matrix.shape[0]):
    for j in range(disease_mobidb_matrix.shape[1]):
        # if disease_mobidb_matrix[i, j] != 0:
        k = int(round(disease_mobidb_matrix[i, j] * 10))
        matrix[i, j, k] = 1

## Dictionary Homo sapiens
for each_feature in mobidb_features_lst:
    cont_fra_temp_lst = mobidb_transposed_df[each_feature].tolist()
    mobidb_predictors_cont_fra_dict[each_feature] = cont_fra_temp_lst

## Plot for homosapiens
for each_feature in mobidb_features_lst[
                    1:]:  # cuz 1st item is 'acc' and we don't need it for the plots, just need the content fraction
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
