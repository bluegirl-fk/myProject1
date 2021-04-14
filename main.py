#Start date = April 14th
# Try to use pivot_table to process the input TSV as shown in the page below
# https://stackoverflow.com/questions/46375147/create-new-columns-from-unique-row-values-in-a-pandas-dataframe
# df = df.pivot_table(index=['acc'], columns=['feature'], values='content_fraction').fillna(0)
#TODO:
import pandas as pd
#TODO: explain each variable and abbreviations in the beginning, define your standard naming system
#TODO: methods
#TODO: import your files
mobidb_original_df = pd.read_csv('data/mobidb_result.tsv', sep='\t')
#transpose the mobidb_original_df and keep only the necessary data
mobidb_transposed_df = mobidb_original_df.pivot_table(index=['acc'], columns=['feature'], values='content_fraction').fillna(0)








#TODO: filter you file with different feature
#TODO: draw each plot
#TODO: Draw comparative plots
#TODO: other stuff
