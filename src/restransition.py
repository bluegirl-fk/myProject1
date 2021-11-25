import pandas as pd
import config as cfg
import seaborn as sns
import matplotlib.pyplot as plt


def residue_heatmapper(df_lst, hmap_title_lst, filename):
    # input df is disorder_majority or filtered_dis_maj if preferred
    new_pivot_df_lst, aa_symbols_lst = heatmap_pivotdf_maker(df_lst)
    sns.set(font_scale=2.7)
    fig, axes = plt.subplots(len(new_pivot_df_lst), 1, figsize=(30, 20 * len(new_pivot_df_lst)))
    for i, (ax, d, t) in enumerate(zip(axes.reshape(-1), new_pivot_df_lst, hmap_title_lst)):
        sb = sns.heatmap(d, cmap="viridis", annot=True, fmt='g', linewidth=0.9, ax=ax, square=True,
                         cbar_kws={'label': 'Transition percentage'}
                         )
        ax.set_title(t, fontsize=40)
        ax.set_xlabel('Variant residues', fontsize=30)
        ax.set_ylabel('Original residues', fontsize=30)
        ax.tick_params(axis='x', colors='red')
        ax.set_xticklabels(aa_symbols_lst)

        # if i < (len(new_pivot_df_lst) - 1):
        #     sb.set(xticklabels=[])
        #     sb.set(xlabel=None)
    plt.tight_layout()
    plt.savefig(cfg.plots['var-hms'] + '/' + filename + '.png', dpi=120)
    plt.show()

    # ax.set_title(hmap_title, fontsize=40)
    # plt.xlabel('Variant residues', fontsize=25)
    # # ax.xaxis.label.set_color('red')
    # plt.ylabel('Original residues', fontsize=25)
    # # ax.yaxis.label.set_color('blue')
    return


def heatmap_pivotdf_maker(df_lst):
    # todo check if the transitions are counted correctly, because the numbers here are different
    new_pivot_df_lst = []
    for df in df_lst:
        res_df = df[['acc', 'var_id', 'orig_aa', 'var_aa']]
        aa_lst = res_df['orig_aa'].unique().tolist()
        res_df = res_df.loc[res_df.var_aa.isin(aa_lst)]
        residue_ser = res_df.groupby(['orig_aa', 'var_aa']).size().to_frame(name='size').reset_index()
        residue_pivot_df = pd.pivot_table(residue_ser, values='size', index='orig_aa', columns='var_aa')
        aa_categ_order_x = ['C', 'M', 'V', 'L', 'I', 'W', 'Y', 'F', 'H', 'T', 'N', 'Q', 'K', 'R', 'D', 'E', 'A', 'S',
                            'G',
                            'P']
        aa_categ_order_y = aa_categ_order_x.__reversed__()
        residue_pivot_df = residue_pivot_df.reindex(columns=aa_categ_order_x)
        residue_pivot_df = residue_pivot_df.reindex(aa_categ_order_y)
        max_value = residue_ser['size'].max()
        print(max_value)
        residue_pivot_df = (residue_pivot_df.div(max_value)).round(3) * 100
        residue_pivot_df = residue_pivot_df.apply(pd.to_numeric)
        new_pivot_df_lst.append(residue_pivot_df)
    difference_df = new_pivot_df_lst[0].subtract(new_pivot_df_lst[1])
    new_pivot_df_lst.append(difference_df)
    return new_pivot_df_lst, aa_categ_order_x


if __name__ == '__main__':
    ## NDD protein list
    ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
    ndd_subdf = ndd_subdf.drop_duplicates()  # (4531, 3)
    ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins
    ## Disorder majority
    # disorder_majority = pd.read_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv')
    # dm_vars_in = disorder_majority.loc[disorder_majority['isin_idr'] == 1]  # (8120, 18)
    # dm_vars_out = disorder_majority.loc[disorder_majority['isin_idr'] == 0]  # (37672, 18)
    # ndd_dm_vars_in = dm_vars_in.loc[dm_vars_in.acc.isin(ndd_pr_lst)]
    # ndd_dm_vars_out = dm_vars_out.loc[dm_vars_out.acc.isin(ndd_pr_lst)]
    ## mobidb_lite
    # 4631 unique ACCs
    mobidb_lite = pd.read_csv(cfg.data['vars'] + '/mobidb_lite-inout-idr-vars-count-normalized.csv')
    mobilite_vars_in = mobidb_lite.loc[mobidb_lite['isin_idr'] == 1]  # (8120, 18)
    mobilite_vars_out = mobidb_lite.loc[mobidb_lite['isin_idr'] == 0]  # (37672, 18)
    ndd_mobilite_vars_in = mobilite_vars_in.loc[mobilite_vars_in.acc.isin(ndd_pr_lst)]
    ndd_mobilite_vars_out = mobilite_vars_out.loc[mobilite_vars_out.acc.isin(ndd_pr_lst)]

    ## heatmaps
    ndd_mobilite_vars_in = ndd_mobilite_vars_in.reset_index()
    res_df = ndd_mobilite_vars_in[['acc', 'var_id', 'orig_aa', 'var_aa']]
    aa_categ_order_x = ['C', 'M', 'V', 'L', 'I', 'W', 'Y', 'F', 'H', 'T', 'N', 'Q', 'K', 'R', 'D', 'E', 'A', 'S', 'G',
                        'P']
    res_df = res_df.loc[res_df.var_aa.isin(aa_categ_order_x)]  # this keeps point transitions and dels the aa additions
    residue_ser = res_df.groupby(['orig_aa', 'var_aa']).size().to_frame(name='size').reset_index()
    residue_pivot_df = pd.pivot_table(residue_ser, values='size', index='orig_aa', columns='var_aa')

    aa_categ_order_y = aa_categ_order_x.__reversed__()
    residue_pivot_df = residue_pivot_df.reindex(columns=aa_categ_order_x)
    residue_pivot_df = residue_pivot_df.reindex(aa_categ_order_y)
    max_value = residue_ser['size'].max()
    print(max_value)
    residue_pivot_df = (residue_pivot_df.div(max_value)).round(3) * 100
    residue_pivot_df = residue_pivot_df.apply(pd.to_numeric)

    # todo turn this into a method to be used in varanalysis