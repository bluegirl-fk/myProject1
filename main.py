# Start date = April 14th
import itertools

import matplotlib.pyplot
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import sys

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
    g.set_xticklabels(xticklabel, rotation=45, va="center", position=(0, -0.02))
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


def matrix_maker_nan(input_df, max_value, thrd_dim_cells, math_oper, get_values_in_range):
    # gets input_df and gives 2d, 3d arrays, filters based on max_value, the 3rd dim is based on the range of values,
    # e.g: 0-10 (thrd_dim_cells) math_oper and get_values_in_range are dependant, they should be calculated in order
    # to give distribute values 0-max based on number of 3d_cells example, for length, max_v = 1000 and we have 11
    # cell ranges [0:11], so we should divide all by 100 => 1000/100 = 10)
    matrix_2d = (input_df.to_numpy() <= max_value) * input_df.to_numpy()
    matrix_3d = np.full((matrix_2d.shape[0], matrix_2d.shape[1], thrd_dim_cells), np.nan)
    for i in range(matrix_2d.shape[0]):
        for j in range(matrix_2d.shape[1]):
            if matrix_2d[i, j] != 0:
                if math_oper == '*':
                    k = int(round(matrix_2d[i, j] * get_values_in_range))
                    matrix_3d[i, j, k] = 1
                elif math_oper == '/':
                    k = int(round(matrix_2d[i, j] / get_values_in_range))
                    matrix_3d[i, j, k] = 1

    #
    # Replace NaN with zeros for rows containing at least one value
    # for i in range(matrix_3d.shape[0]):
    #     for j in range(matrix_3d.shape[1]):
    #         if 1.0 in matrix_3d[i][j]:
    #             matrix_3d[i][j][np.isnan(matrix_3d[i][j])] = 0

    # sum of Pr.s with same content_fraction or length for each feature
    matrix_3d_sum = np.nansum(matrix_3d, axis=0, dtype=int)
    matrix_3d_sum_normalized = matrix_3d_sum / matrix_3d_sum.max(axis=1)[:, None]
    return matrix_2d, matrix_3d, matrix_3d_sum, matrix_3d_sum_normalized


def matrix_maker_zeros(input_df, num_3rd_dim, max_range):
    matrix_2d = (input_df.to_numpy() <= 1.) * input_df.to_numpy()
    matrix_3d = np.zeros((matrix_2d.shape[0], matrix_2d.shape[1], num_3rd_dim + 1))
    for i in range(matrix_2d.shape[0]):
        for j in range(matrix_2d.shape[1]):
            # if matrix_2d[i, j] != 0:
            k = int(round(matrix_2d[i, j] * num_3rd_dim))
            matrix_3d[i, j, k] = 1
    matrix_3d_sum = np.sum(matrix_3d, axis=0)
    # matrix_3d_sum_normalized = matrix_3d_sum / matrix_3d_sum.max(axis=1)[:, None]
    return matrix_2d, matrix_3d, matrix_3d_sum  # , matrix_3d_sum_normalized


def sum_df_generator(input_sum_matrix, columns):
    sum_df = pd.DataFrame(input_sum_matrix,
                          columns=columns,
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


def df_lst_maker_for_barplot(input_matrix):
    cols_sum_df = pd.DataFrame([input_matrix], columns=mobidb_features_lst[1:], index=['Protein count'])
    cols_sum_df = cols_sum_df.T.reset_index()
    cols_sum_df.columns = ['Features', 'Protein count']
    cols_sum_lst = cols_sum_df['Protein count']
    cols_sum_lst = [int(x) for x in cols_sum_lst]

    return cols_sum_df, cols_sum_lst


def expand_regions(region_ranges_lst):
    transformed_regions = []
    for reg in region_ranges_lst:
        start = int(reg.split('..')[0])
        end = int(reg.split('..')[1])
        while start <= end:
            transformed_regions.append(start)
            start += 1
    return set(transformed_regions)


if __name__ == '__main__':
    # This is the main!

    ## Gene4denovo
    # g4dn with only exonic mutations
    exonic_g4dn_df = pd.read_csv('data/gene4denovo/exonic-df.csv')  # (70879, 156)
    del exonic_g4dn_df['Unnamed: 0']
    # set its own index to idx, npot the original df , then m,erge with subdf
    exonic_g4dn_df = exonic_g4dn_df.reset_index()
    exonic_g4dn_df = exonic_g4dn_df.rename(columns={'index': 'idx', 'AAChange.refGene': 'AAChange_refGene'})

    # aachange_g4dn_subdf1 = exonic_g4dn_df[['AAChange_refGene']]
    # aachange_g4dn_subdf1 = aachange_g4dn_subdf1['AAChange_refGene'].str.split(',', expand=True).stack().to_frame(
    #     'AAChange_refGene')
    # aachange_g4dn_subdf2 = aachange_g4dn_subdf1['AAChange_refGene'].str.split(pat=':', expand=True)
    # aachange_g4dn_subdf3 = aachange_g4dn_subdf2[aachange_g4dn_subdf2.columns[:10]]  # del rest of None columns
    # aachange_g4dn_subdf3.columns = ['Gene_refGene', 'refSeq', 'exon#', 'mutNA', 'AAChange_refGene', 'aa1', 'aa2', 'position',
    #                          'frameshift', 'mutPr']  # rename cuz .stack() made Series and deleted original column names
    # aachange_g4dn_subdf3['mutPr'] = aachange_g4dn_subdf3.AAChange_refGene.str.split(pat='fs*', expand=True, )
    # # sys.exit()
    # aachange_g4dn_subdf3['aa1'] = aachange_g4dn_subdf3['mutPr'].str[2]
    # aachange_g4dn_subdf3['aa2'] = aachange_g4dn_subdf3['mutPr'].str[-1]
    # aachange_g4dn_subdf3['position'] = aachange_g4dn_subdf3['mutPr'].str.replace(r'\D', '')  # (201372, 10)
    # aachange_g4dn_subdf3['frameshift'] = aachange_g4dn_subdf3['AAChange_refGene'].str.split('fs', 1).str[1]  # find out
    # # scientific meaning of fs*35 for example
    # del aachange_g4dn_subdf3['mutPr']
    # aachange_g4dn_subdf3.to_csv(r'data/gene4denovo/subdf-mut-beforeACC.csv')  # (201372, 9)
    # # Then made a list of refSeq ids from this df, wrote to .txt, retrived ACCs from uniprot
    # # splited my text file using bash : split -l 70000 refseq-gene4dn.txt, the 7000 is number of the lines
    # TODO: filter phens

    ## mut positionb df import
    rseq_mutinfo_df = pd.read_csv('data/gene4denovo/subdf-mut-beforeACC.csv')  # (201372, 11)
    rseq_mutinfo_df = rseq_mutinfo_df.rename(columns={'Unnamed: 0': 'idx', 'Unnamed: 1': 'sub-idx'})
    # merge with gene4dn exonic, now original exonic g4dn file + mutInfo e.g position
    g4dn_exonic_mutinfo_df = pd.merge(rseq_mutinfo_df, exonic_g4dn_df, on='idx')  # (201372, 166)
    ## g4dn mutInfo + uniprot ACCs file
    refseq_acc_df1 = pd.read_csv('data/refseq/refseq-acc.tab', sep='\t')  # from Uniprot
    refseq_acc_df1.columns = ['refseq_id', 'isoforms', 'acc', 'organism', 'Length', 'Gene names']
    del refseq_acc_df1['organism']  # (50930, 5)
    refseq_acc_df2 = refseq_acc_df1['refseq_id'].str.split(',', expand=True).stack()
    idx_tmp = refseq_acc_df2.index._get_level_values(0)
    refseq_acc_df3 = refseq_acc_df1.iloc[idx_tmp].copy()
    refseq_acc_df3['refSeq'] = refseq_acc_df2.values
    del refseq_acc_df3['refseq_id']  # (93949, 5)

    # merge with aachange_g4dn_subdf3 to get positions
    # rseq_mutinfo_df = pd.read_csv('data/gene4denovo/subdf-mut-beforeACC.csv')  # (201372, 11)
    # rseq_mutinfo_df = rseq_mutinfo_df.rename(columns={'Unnamed: 0': 'idx', 'Unnamed: 1': 'sub-idx'})
    mut_positions_df = pd.merge(refseq_acc_df3, rseq_mutinfo_df, on='refSeq')
    del mut_positions_df['Gene names']
    mut_positions_df = mut_positions_df.loc[:, ~mut_positions_df.columns.str.contains('^Unnamed')]
    mut_positions_df = mut_positions_df[
        ['acc', 'position', 'aa1', 'aa2', 'Gene.refGene', 'refSeq', 'Length', 'isoforms', 'exon#', 'mutNA', 'AAChange_refGene',
         'frameshift']]
    mut_positions_df['position'] = mut_positions_df['position'].fillna(0).astype(int)  # (551773, 12)
    mut_positions_df = mut_positions_df.drop_duplicates(ignore_index=True)  # (224021, 12)
    mut_positions_df.reset_index(level=0, inplace=True)

    ## mobidb
    mobidb_original_df = pd.read_csv('data/mobidb_result.tsv', sep='\t')
    mobidb_original_df.columns = ['acc', 'feature', 'startend', 'content_fraction', 'content_count', 'length']
    # converting disorder content ranges in each cell to a list
    mobidb_original_df['startend'] = mobidb_original_df['startend'].str.split(',')
    mut_pos_subdf = mut_positions_df[['index', 'acc', 'position']]

    final_mut_check_df = pd.read_csv('data/mutations-position-mobidb.csv')
    filtered_mut_pos_df = final_mut_check_df[final_mut_check_df['is_in_startend'] == 1]  # (941013, 10)
    unique_mut_pos_df = filtered_mut_pos_df.drop_duplicates(['acc', 'startend', 'position', 'length'])
    # Merge the dataframes
    mobidb_mutpos_df = pd.merge(mobidb_original_df, mut_pos_subdf, on='acc')  # (4013158, 8)
    array_is_in = []
    for index, row in mobidb_mutpos_df.iterrows():
        set_disorder_region = expand_regions(row.startend)
        print(set_disorder_region)
        if row.position in set_disorder_region:
            # print(mobidb_mutpos_df.loc[mobidb_mutpos_df['position'], 'index'])
            print('yes')
            array_is_in.append('1')
        else:
            print('no')
            array_is_in.append('0')

    mobidb_mutpos_df['is_in_startend'] = array_is_in
    mobidb_mutpos_df.to_csv(r'data/mutations-position-mobidb.csv')
    # TODO: filter the phenotypes from gene4dn, like the phens file


    ## DBsnp
    # with uniprot acc
    snpdb_acc = pd.read_csv('data/refsnp/uniprot.tsv', sep='\t')  # (10219, 3)
    snpdb_acc.columns = ['avsnp150', 'rs_ids', 'uniprot_acc']
    snpdb_acc[['acc', 'variant']] = snpdb_acc.uniprot_acc.str.split('#', expand=True)
    snpdb_acc = snpdb_acc.drop(columns=['rs_ids', 'uniprot_acc'])
    # snpdb_merged = pd.merge(snpdb_acc, snpdb, on='avsnp150') # (158780,8) # problem?
    # for now don't use the merged df, uniprot not for now, ensembel alternative

    snpdb = pd.read_csv('data/refsnp/genebank.tsv', sep='\t')  # (946889, 6)
    # snpdb.columns = ['avsnp150', 'rs_ids', 'seq_id', 'position', 'del_seq', 'in_seq']  # avsnp150 = refsnp_id
    # snpdb_mut = snpdb[snpdb.del_seq != snpdb.in_seq]  # (338184, 6)
    # snpdb_mut = snpdb_mut.drop_duplicates(ignore_index=True)  # (327147, 7)
    # snpdb_mut['all_rsids'] = snpdb_mut[['avsnp150', 'rs_ids']].astype(str).agg(','.join, axis=1)
    # snpdb_mut.to_csv(r'data/snpdb_mut.tsv')

    snpdb_mut = pd.read_csv('data/dbsnp_mutations.tsv')
    snpdb_mut_dict = dict(zip(snpdb_mut.index, snpdb_mut.all_rsids))  # 327147
    keys_values_snpdb_mut_dict = snpdb_mut_dict.items()
    snpdb_mut_dict = {str(key): list(str(value).split(",")) for key, value in keys_values_snpdb_mut_dict}  # dict with
    # idx of snpdb_mut df as key and lst of all merged rsids of each position as values

    snpdb_idx_true_lst = []
    snpdb_idx_false_lst = []  # not used, can come handy for not mapped ones if needed

    # this for loop checks if rs_ids of g4dn are in our filtered snpdb dataset or not. then idxs are added to .txt file
    # for i in g4dn_rsid_lst:  # is this correct? why less numbers than merging two dfs?
    #     for key, values in snpdb_mut_dict.items():
    #         if (isinstance(values, list)):
    #             if i in values:
    #                 snpdb_idx_true_lst.append(key)  # len : 984

    true_snpdb_g4dn_idxs = pd.read_csv('data/dbsnp-idxs.txt')
    snpdb_idx_true_lst = true_snpdb_g4dn_idxs['index'].tolist()

    snpdb_idx_true_df = snpdb_mut.loc[snpdb_mut.index[snpdb_idx_true_lst]]  # (984, 7), 57 unique rsids
    # uniprot_db_g4dn = pd.merge(snpdb_acc, snpdb_g4dn_df, on='avsnp150')  # only contains 1/2

    snpdb_idx_true_df['avsnp150'] = snpdb_idx_true_df['avsnp150'].astype(int)
    g4dn_df['avsnp150'] = g4dn_df['avsnp150'].astype(int)
    snpdb_mut['avsnp150'] = snpdb_mut['avsnp150'].astype(int)

    g4dn_mapped_snpdb_df = pd.merge(snpdb_idx_true_df, g4dn_df,
                                    on='avsnp150')  # (1000, 13) # trying to get snpdb_idx_true_df + columns
    # of gene4dn e.g: esembel id, etc
    g4dn_snpdb_acc_df = pd.merge(g4dn_mapped_snpdb_df, snpdb_acc, on='avsnp150')  # (410, 15)

    snpdb_acc_pos_df = g4dn_snpdb_acc_df.drop_duplicates(['acc', 'position'])  # (26, 15)

    snpdb_g4dn_idx_merge_df = pd.merge(g4dn_df, snpdb_mut, on='avsnp150')  # (976, 13)
    # positions in my file should be +1 for the equivalent one in uniprot

    ### Files import and modify

    ## for content fraction
    mobidb_pivot_contf_df = mobidb_original_df.pivot_table(
        index=['acc'], columns=['feature'], values='content_fraction').fillna(0)
    mobidb_pivot_contf_df = mobidb_pivot_contf_df.reset_index()  # added idx nums manually,ACCs recognized separate column
    mobidb_pivot_contf_df.to_csv(r'data/mobidb_pivot_contf_df.csv', index=True)
    mobidb_features_lst = mobidb_pivot_contf_df.columns.str.split(',').tolist()  # this also contains the 'acc' column
    mobidb_features_lst = list(itertools.chain(*mobidb_features_lst))  # flat list

    ndd_acc_df = pd.read_csv('data/allUniqueEntry.tab', sep='\t')
    ndd_acc_lst = ndd_acc_df['Entry'].to_list()
    ndd_contf_df = mobidb_pivot_contf_df[mobidb_pivot_contf_df['acc'].isin(ndd_acc_lst)]
    ndd_contf_df.to_csv(r'data/ndd-contf-dataframe.csv')
    ## for Length
    mobidb_length_df = mobidb_original_df[['acc', 'length']].drop_duplicates(subset=['acc'])
    mobidb_pivot_length_df = mobidb_original_df.pivot_table(
        index=['acc'], columns=['feature'], values='length').fillna(0)
    mobidb_pivot_length_df = mobidb_pivot_length_df.reset_index()  # reset idx to get acc as dif col to search in it
    mobidb_pivot_length_df.to_csv(r'data/mobidb-pivot-length-df.csv')
    ndd_length_df = mobidb_pivot_length_df[mobidb_pivot_length_df['acc'].isin(ndd_acc_lst)]  # len(df) = 1089
    ndd_length_df.to_csv(r'data/ndd-length-df.csv')
    ## Matrix
    # content fraction with nan
    _, mobi_contf_mat, mobi_contf_sum_mat, mobi_contf_sum_norm_mat = matrix_maker_nan(
        input_df=mobidb_pivot_contf_df.iloc[:, 1:], max_value=1., thrd_dim_cells=11, math_oper='*',
        get_values_in_range=10)
    _, ndd_contf_mat, ndd_contf_sum_mat, ndd_contf_sum_norm_mat = matrix_maker_nan(
        input_df=ndd_contf_df.iloc[:, 1:], max_value=1., thrd_dim_cells=11, math_oper='*', get_values_in_range=10)

    # Length (Use vstack or hstack ?)
    mobi_len_2d_mat, mobi_len_3d_mat, mobi_len_sum_mat, mobi_len_sum_norm_mat = matrix_maker_nan(
        input_df=mobidb_pivot_length_df.iloc[:, 1:], max_value=1000, thrd_dim_cells=11, math_oper='/',
        get_values_in_range=100)
    ndd_len_2d_mat, ndd_len_3d_mat, ndd_len_sum_mat, ndd_len_sum_norm_mat = matrix_maker_nan(
        input_df=ndd_length_df.iloc[:, 1:], max_value=1000, thrd_dim_cells=11, math_oper='/', get_values_in_range=100)
    # ax.hist(dataset_len, bins=np.arange(0, 1000, 10))

    ## sum dataframes
    mobi_contf_sum_norm_df = sum_df_generator(mobi_contf_sum_norm_mat,
                                              ['0', ' ', '20', ' ', '40', ' ', '60', ' ', '80', ' ', '100'])
    ndd_cont_fract_sum_norm_df = sum_df_generator(ndd_contf_sum_norm_mat,
                                                  ['0', ' ', '20', ' ', '40', ' ', '60', ' ', '80', ' ', '100'])
    mobi_len_sum_norm_df = sum_df_generator(mobi_len_sum_norm_mat,
                                            [' ', '100', ' ', '300', ' ', '500', ' ', '700', ' ', '900', ''])
    ndd_len_sum_norm_df = sum_df_generator(ndd_len_sum_norm_mat,
                                           [' ', '100', ' ', '300', ' ', '500', ' ', '700', ' ', '900', ''])

    ## Difference of the sum arrays(with nan)
    difference_contf_sum_norm_mat = mobi_contf_sum_norm_mat - ndd_contf_sum_norm_mat
    difference_len_sum_norm_mat = mobi_len_sum_norm_mat - ndd_len_sum_norm_mat
    difference_contf_sum_norm_df = sum_df_generator(difference_contf_sum_norm_mat,
                                                    ['0', ' ', '20', ' ', '40', ' ', '60', ' ', '80', ' ', '100'])
    difference_len_sum_norm_df = sum_df_generator(difference_len_sum_norm_mat,
                                                  [' ', '100', ' ', '300', ' ', '500', ' ', '700', ' ', '900', ''])

    ## heatmaps
    draw_heatmaps([mobi_contf_sum_norm_df.T, ndd_cont_fract_sum_norm_df.T, difference_contf_sum_norm_df.T],
                  ['Homo sapiens', 'NDDs', 'Difference (Homo sapiens - NDDs)'],
                  saving_rout='plots/heatmaps/Heatmaps0.png')
    draw_heatmaps([mobi_len_sum_norm_df.T, ndd_len_sum_norm_df.T, difference_len_sum_norm_df.T],
                  ['Homo sapiens', 'NDDs', 'Difference (Homo sapiens - NDDs)'],
                  saving_rout='plots/heatmaps/Heatmaps-Length.png')

    ## columns sum sum_matrix to get protein numbers (for bar plot based on distribution of heatmap)
    # content fraction
    mobi_contf_cols_sum_df, mobi_contf_cols_sum_lst = df_lst_maker_for_barplot(mobi_contf_sum_mat.T.sum(axis=0))
    ndd_contf_cols_sum_df, ndd_contf_cols_sum_lst = df_lst_maker_for_barplot(ndd_contf_sum_mat.T.sum(axis=0))
    mobi_contf_cols_sum_df.to_csv(r'data/mobidb-contf-distribution.csv')
    ndd_contf_cols_sum_df.to_csv(r'data/ndd-contf-distribution.csv')
    # length
    mobi_len_cols_sum_df, mobi_len_cols_sum_lst = df_lst_maker_for_barplot(mobi_len_sum_mat.T.sum(axis=0))
    ndd_len_cols_sum_df, ndd_len_cols_sum_lst = df_lst_maker_for_barplot(ndd_len_sum_mat.T.sum(axis=0))
    mobi_len_cols_sum_df.to_csv(r'data/mobidb-len-distribution.csv')
    ndd_len_cols_sum_df.to_csv(r'data/ndd-len-distibution.csv')

    ## Hmap distribution barplots (Protein count)
    # content fraction
    draw_barplot(figsize_a='18', figsize_b='9', xlabel='Features', ylabel='Protein count', data=mobi_contf_cols_sum_df,
                 xticklabel=mobi_contf_cols_sum_lst, yscale='log',
                 save_rout='plots/log/hist-hmaps-distribution/mobidb-contf-NEW.png')
    draw_barplot(figsize_a='18', figsize_b='9', xlabel='Features', ylabel='Protein count', data=ndd_contf_cols_sum_df,
                 xticklabel=ndd_contf_cols_sum_lst, yscale='log',
                 save_rout='plots/log/hist-hmaps-distribution/ndd-contf-NEW.png')
    # Length
    draw_barplot(figsize_a='18', figsize_b='9', xlabel='Features', ylabel='Protein count', data=mobi_len_cols_sum_df,
                 xticklabel=mobi_len_cols_sum_lst, yscale='log',
                 save_rout='plots/log/hist-hmaps-distribution/mobi-len-NEW.png')
    draw_barplot(figsize_a='18', figsize_b='9', xlabel='Features', ylabel='Protein count', data=ndd_len_cols_sum_df,
                 xticklabel=ndd_len_cols_sum_lst, yscale='log',
                 save_rout='plots/log/hist-hmaps-distribution/ndd-len-NEW.png')

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
        ndd_cont_fra_temp_lst = ndd_contf_df[each_feature].tolist()
        ndd_predictors_cont_fra_dict[each_feature] = ndd_cont_fra_temp_lst

    ## plot for ndds
    for each_feature in mobidb_features_lst[1:]:
        drawplot(ndd_predictors_cont_fra_dict[each_feature], 'log', 30, False, each_feature + '_NDD', 'Protein count',
                 'plots/log/hist-all-NDD/SCZ/' + each_feature + '_SCZ')

    # comparative histogram (homosapiens Vs. ndd)
    for each_feature in mobidb_features_lst[1:]:
        compare_plot(first_lst=mobidb_predictors_cont_fra_dict[each_feature],
                     second_lst=ndd_predictors_cont_fra_dict[each_feature], yscale='log', bins=30, is_dense=False,
                     x_label=each_feature + '_comparison', y_label='proteins count',
                     first_label='Homo sapiens Pr.s', second_label='NDD Pr.s',
                     png_file_name='plots/log/hist-comparison' '-homoS-NDD/' + each_feature)
