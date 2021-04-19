# Start date = April 14th
### import
import itertools

import pandas as pd
import matplotlib.pyplot as plt

### explain variables and stuff
mobidb_features_lst = []  # setting the transversed dataframe's column names which are the feature/predictor names as
# list values,#column names list contains dif features (dif predictors for dif plots, 78 features = 78 plots)
mobidb_predictors_cont_fra_dict = {}  # dict with each feature as key and the associated list of content_fractions as
# values
cont_fra_temp_lst = []  # # content_fraction of homo sapiens list based on each feature,


# it changes with the for loop iteration


### methods
def drawplot(plot_input_lst, bins_amount, isDense, xLabel, yLabel, png_file_name):
    plt.hist(plot_input_lst, bins=bins_amount,
             density=isDense)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    file_format = '.png'
    plt.savefig('graphs/' + png_file_name + '_content_fraction' + file_format)
    plt.show()



### import files
mobidb_original_df = pd.read_csv('data/mobidb_result.tsv', sep='\t')
# transpose the mobidb_original_df and keep only the necessary data
mobidb_transposed_df = mobidb_original_df.pivot_table(index=['acc'], columns=['feature'],
                                                      values='content_fraction').fillna(0)
mobidb_transposed_df.to_csv(r'data/mobidb_transposed_df.csv', index=True)
### extract required data from the dataframe
mobidb_features_lst = mobidb_transposed_df.columns.str.split(',').tolist()
mobidb_features_lst = list(itertools.chain(*mobidb_features_lst))  # list of lists to a flat list
# get each column's content as a list and set it as dict value (keys are predictor's name(feature's name))
for each_feature in mobidb_features_lst:
    cont_fra_temp_lst = mobidb_transposed_df[each_feature].tolist()
    mobidb_predictors_cont_fra_dict[each_feature] = cont_fra_temp_lst
### draw each plot
for each_feature_name in mobidb_features_lst:
    drawplot(mobidb_predictors_cont_fra_dict[each_feature_name], 10, True, each_feature_name, "Protein Count",
             each_feature_name)  # maybe need to use str() for the file name and x axis label

# TODO: Draw comparative plots

# TODO: check if it is needed to add zeros and sum = 75088 instead of 75052 //or omitted data:/ too much zeros
