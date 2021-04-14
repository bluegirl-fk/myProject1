#Start date = April 14th
### import
import pandas as pd
import matplotlib.pyplot as plt

#TODO: explain each variable and abbreviations in the beginning, define your standard naming system

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
### extract required data from the datafram
mobidb_transposed_df_features_lst = mobidb_transposed_df.columns.str.split(',').tolist() #column names list contains dif features (dif predictors for dif plots, 78 features = 78 plots)

#TODO: filter you file with different feature
#find content_fraction for each of the features iterating through all the proteins list


#TODO: draw each plot
#TODO: Draw comparative plots
#TODO: other stuff

