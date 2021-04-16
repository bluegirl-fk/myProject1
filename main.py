#Start date = April 14th
### import
import itertools

import pandas as pd
import matplotlib.pyplot as plt

#TODO: explain each variable and abbreviations in the beginning, define your standard naming system
mobidb_features_lst = [] #setting the transversed dataframe's column names which are the feature/predictor names as list values
mobidb_predictors_cont_fra_dict = {}
### methods
def drawplot(plot_input_lst, bins_amount, isDense, xLabel, yLabel, png_file_name):
    plt.hist(plot_input_lst, bins=bins_amount,
             density=isDense)  # density=True shows comparison percentage in y axis (that value compared to whole data)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.show()
    plt.savefig(png_file_name.png)

### import your files
mobidb_original_df = pd.read_csv('data/mobidb_result.tsv', sep='\t')
#transpose the mobidb_original_df and keep only the necessary data
mobidb_transposed_df = mobidb_original_df.pivot_table(index=['acc'], columns=['feature'], values='content_fraction').fillna(0)
mobidb_transposed_df.to_csv(r'data/mobidb_transposed_df.csv', index=True)
### extract required data from the dataframe
mobidb_features_lst = mobidb_transposed_df.columns.str.split(',').tolist() #column names list contains dif features (dif predictors for dif plots, 78 features = 78 plots)
mobidb_features_lst = list(itertools.chain(*mobidb_features_lst)) #list of lists to a flat list
#get each column's content as a list and set it as dict value (keys are predictor's name(feature's name))
print(mobidb_features_lst)
for each_feature in mobidb_features_lst:
    mobidb_each_cont_fra_lst = mobidb_transposed_df[each_feature].tolist()  # content_fraction of homo sapiens list based on each feature
    mobidb_predictors_cont_fra_dict[each_feature] = mobidb_each_cont_fra_lst
#TODO: draw each plot
#TODO: Draw comparative plots
#TODO: other stuff

