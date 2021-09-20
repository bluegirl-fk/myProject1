import os

absolute = os.path.abspath(os.getcwd())
absolute = absolute + '/../'
data = {
    'gene4': absolute + 'data/gene4denovo',
    'rseq': absolute + 'data/refseq',
    '': absolute + 'data',
    'brain': absolute + 'data/brain',
    'csv': absolute + 'data/to-csv',
    'dndb': absolute + 'data/denovodb',
    'phens': absolute + 'data/mobibool-phens',
    'phens-fdr': absolute + 'data/mobibool-phens/acc-phens-FDR',
    'xml': absolute + 'data/xml',
    'xml-p': absolute + 'data/xml/parsed',
    'vars': absolute + 'data/variants',
    'var-desc': absolute + 'data/variants/descriptions',
    'var-desc-cf': absolute + 'data/variants/descriptions/cf',
    'var-desc-cc': absolute + 'data/variants/descriptions/cc',
    'var-desc-len': absolute + 'data/variants/descriptions/length'

}
plots = {
    '': absolute + 'plots',
    'bar': absolute + 'plots/g4dn/barplot',
    'box': absolute + 'plots/g4dn/boxplot',
    'box-cf': absolute + 'plots/g4dn/boxplot/cf',
    'box-cc': absolute + 'plots/g4dn/boxplot/cc',
    'box-len': absolute + 'plots/g4dn/boxplot/length',
    'violin': absolute + 'plots/g4dn/violinplot',
    'vio-cf': absolute + 'plots/g4dn/violinplot/cf',
    'vio-cc': absolute + 'plots/g4dn/violinplot/cc',
    'vio-len': absolute + 'plots/g4dn/violinplot/length',
    'vars': absolute + 'plots/variants',
    'vb-cf': absolute + 'plots/variants/boxplot/cf',
    'vb-cc': absolute + 'plots/variants/boxplot/cc',
    'vb-len': absolute + 'plots/variants/boxplot/length',
    'vv-cf': absolute + 'plots/variants/violinplot/cf',
    'vv-cc': absolute + 'plots/variants/violinplot/cc',
    'vv-len': absolute + 'plots/variants/violinplot/length',
    'var-bar': absolute + 'plots/variants/barplot',


}
