import config as cfg
import pandas as pd
import sys

def main():
    generatefff()
    # generateggg()

def generatefff():

    # ### G4dn code more ## Gene4denovo ## only exonic mutations # exonic_g4dn_df = pd.read_csv(
    # 'data/gene4denovo/exonic-df.csv')  # (70879, 156) # stacked the refseq mut positions, now have repeated
    # ...............

    df = pd.read_csv(cfg.data['gene4'] + '/yourfile.csv', sep='\t')
    # ...............

def generateggg():
    # ### G4dn code more ## Gene4denovo ## only exonic mutations # exonic_g4dn_df = pd.read_csv(
    # 'data/gene4denovo/exonic-df.csv')  # (70879, 156) # stacked the refseq mut positions, now have repeated
    # ...............

    df = pd.read_csv(cfg.data['gene4'] + '/yourfile.csv', sep='\t')
    # ...............



if __name__ == '__main__':
    sys.exit(main())