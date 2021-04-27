# Start date = April 14th
### import
import itertools

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

### explain variables and stuff
mobidb_features_lst = []  # setting the transversed dataframe's column names which are the feature/predictor names as
# list values,#column names list contains dif features (dif predictors for dif plots, 78 features = 78 plots)
mobidb_predictors_cont_fra_dict = {}  # dict with each feature as key and the associated list of content_fractions as
# values
cont_fra_temp_lst = []  # content_fraction of homo sapiens list based on each feature, it changes with the for loop iteration

disease_cont_fra_temp_lst = []  # content_fraction of disease acc list based on each feature, it changes with the for loop iteration
disease_predictors_cont_fra_dict = {}  # dict with each feature as key and the associated list of ONLY disease content_fractions as values


### methods
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


if __name__ == '__main__':
    ### Dataframe files import and manipulation

    mobidb_original_df = pd.read_csv('data/mobidb_result.tsv', sep='\t')
    # transpose the mobidb_original_df and keep only the necessary data
    mobidb_transposed_df = mobidb_original_df.pivot_table(index=['acc'], columns=['feature'],
                                                          values='content_fraction').fillna(0)
    mobidb_transposed_df = mobidb_transposed_df.reset_index()  # adding index numbers manually to have the acc as a separate column and not as index
    mobidb_transposed_df.to_csv(r'data/mobidb_transposed_df.csv', index=True)

    disease_acc_df = pd.read_csv('data/diseases.tab', sep='\t')
    # list of protein ACCs from disease acc df
    disease_acc_lst = disease_acc_df['Entry'].to_list()
    # new DF containing only rows with disease protein ACC
    disease_mobidb_df = mobidb_transposed_df[mobidb_transposed_df['acc'].isin(disease_acc_lst)]
    ### extract required data from the mobidb tansposed dataframe
    mobidb_features_lst = mobidb_transposed_df.columns.str.split(',').tolist()  # this also contains the 'acc' column
    disease_mobidb_matrix = (disease_mobidb_df.iloc[:, 1:].to_numpy() <= 1.) * disease_mobidb_df.iloc[:, 1:].to_numpy()
    matrix = np.zeros((disease_mobidb_matrix.shape[0], disease_mobidb_matrix.shape[1], 11))
    for i in range(disease_mobidb_matrix.shape[0]):
        for j in range(disease_mobidb_matrix.shape[1]):
            # if disease_mobidb_matrix[i, j] != 0:
            k = int(round(disease_mobidb_matrix[i, j] * 10))
            matrix[i, j, k] = 1
    pass
    # mobidb_features_lst = list(itertools.chain(*mobidb_features_lst))  # list of lists to a flat list
    # # get each column's content as a list and set it as dict value (keys are predictor's name(feature's name))
    # for each_feature in mobidb_features_lst:
    #     cont_fra_temp_lst = mobidb_transposed_df[each_feature].tolist()
    #     mobidb_predictors_cont_fra_dict[each_feature] = cont_fra_temp_lst
    # ### draw each plot of whole homo sapiens proteins with content fraction based on each feature
    # # for each_feature in mobidb_features_lst[1:]: ##cuz the fist item is 'acc' and we don't need it for the plots, just need the content fraction
    # #     drawplot(mobidb_predictors_cont_fra_dict[each_feature], 30, True, each_feature + '_homosapiens', "Protein Count(relative)",
    # #              each_feature, 'homosapiens')
    # ##plot for diseases
    # for each_feature in mobidb_features_lst:
    #     disease_cont_fra_temp_lst = disease_mobidb_df[each_feature].tolist()
    #     disease_predictors_cont_fra_dict[each_feature] = disease_cont_fra_temp_lst
    #
    # # for each_feature in mobidb_features_lst[1:]:
    # #      drawplot(disease_predictors_cont_fra_dict[each_feature], 30, False, each_feature + '_Disease', 'Protein count',each_feature+'_disease','disease')
    #
    # # ### comparative histogram
    # # for each_feature in mobidb_features_lst[
    # #                     1:]:  # cuz the fist item is 'acc' and we don't need it for the plots, just need the content fraction
    # #     compare_plot(mobidb_predictors_cont_fra_dict[each_feature], disease_predictors_cont_fra_dict[each_feature], 30,
    # #                  True, x_label=each_feature + '_comparison', y_label='proteins count(relative)', png_file_name=each_feature,
    # #                  first_label= 'all_disordered_Pr.s', second_label= 'disease_Pr.s')
    #
