import os

absolute = os.path.abspath(os.getcwd())
absolute = absolute + '/../'
data = {
    'gene4': absolute + 'data/gene4denovo',
    'rseq': absolute + 'data/refseq',
    '': absolute + 'data',
    'brain': absolute + 'data/brain',
    'csv': absolute + 'data/to-csv',

}
plots = {
    'box': absolute + 'plots/g4dn/boxplot',
    'violin': absolute + 'plots/g4dn/violinplot',

}
