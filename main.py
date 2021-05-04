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


# TODO: make matrix maker methods, one with nan as Damiano did, one without like you did before

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

## 3d matrix mobidb
mobidb_matrix = (mobidb_transposed_df.iloc[:, 1:].to_numpy() <= 1.) * mobidb_transposed_df.iloc[:, 1:].to_numpy()
# Replace "np.nan" with 0 to initialize the full matrix with zeros
mobidb_3d_matrix = np.full((mobidb_matrix.shape[0], mobidb_matrix.shape[1], 11), np.nan)
for i in range(mobidb_matrix.shape[0]):
    for j in range(mobidb_matrix.shape[1]):
        if mobidb_matrix[i, j] != 0:
            k = int(round(mobidb_matrix[i, j] * 10))
            mobidb_3d_matrix[i, j, k] = 1

# Replace NaN with zeros for rows containing at least one value
for i in range(mobidb_3d_matrix.shape[0]):
    for j in range(mobidb_3d_matrix.shape[1]):
        if 1.0 in mobidb_3d_matrix[i][j]:
            mobidb_3d_matrix[i][j][np.isnan(mobidb_3d_matrix[i][j])] = 0

# gives us sum of content_fraction of all given proteins based on each mobidb feature
mobidb_3d_matrix_sum = np.nansum(mobidb_3d_matrix, axis=0)
mobid_cont_fract_sum_df = pd.DataFrame(mobidb_3d_matrix_sum,
                                       columns=['0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'],
                                       index=mobidb_features_lst[1:])

## 3d matrix disease
disease_mobidb_matrix = (disease_mobidb_df.iloc[:, 1:].to_numpy() <= 1.) * disease_mobidb_df.iloc[:, 1:].to_numpy()
disease_3d_matrix = np.full((disease_mobidb_matrix.shape[0], disease_mobidb_matrix.shape[1], 11), np.nan)
for i in range(disease_mobidb_matrix.shape[0]):
    for j in range(disease_mobidb_matrix.shape[1]):
        if disease_mobidb_matrix[i, j] != 0:
            k = int(round(disease_mobidb_matrix[i, j] * 10))
            disease_3d_matrix[i, j, k] = 1

for i in range(disease_3d_matrix.shape[0]):
    for j in range(disease_3d_matrix.shape[1]):
        if 1.0 in disease_3d_matrix[i][j]:
            disease_3d_matrix[i][j][np.isnan(disease_3d_matrix[i][j])] = 0

disease_3d_matrix_sum = np.nansum(disease_3d_matrix, axis=0)
disease_3d_matrix_sum = disease_3d_matrix_sum / disease_3d_matrix_sum.max(axis=1)[:,None]
disease_cont_fract_sum_df = pd.DataFrame(disease_3d_matrix_sum,
                                         columns=['0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'],
                                         index=mobidb_features_lst[1:])
## heatmaps of Sum
sns.set()
fig, ax = plt.subplots(2, 1, figsize=(24, 12))
for i, d in enumerate([mobid_cont_fract_sum_df.transpose(), disease_cont_fract_sum_df.transpose()]):
    labels = d.applymap(lambda v: str(v) if v == d.values.max() else '')
    sns.heatmap(d,
                cmap="viridis",  # sequential colormap
                annot_kws={'fontsize': 11},
                fmt='',
                square=True,
                #vmax=1,
                #vmin=0,
                linewidth=0.01,
                linecolor="#222",
                ax=ax[i],
                )

ax[0].set_title('Homo sapiens')
ax[1].set_title('NDDs')
ax[0].set_ylabel('Disorder %')
ax[1].set_ylabel('Disorder %')
ax[0].set_xlabel('mobiDB features')
ax[1].set_xlabel('mobiDB features')

plt.tight_layout()
plt.savefig('plots/heatmaps/final.png', dpi=120)
plt.show()

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
