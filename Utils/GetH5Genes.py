import scanpy as sp
import argparse

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-h5', '--h5file', type=str, help = 'H5 file to harvest genes from')
    args = parser.parse_args()

    return args


#This function is designed to try and handle a more broad range of H5 files
#We aim to find where expression data and gene names are stored.
def prep_h5(infile):

    adata = sp.read_h5ad(infile)

    frame_setting = adata

    if adata.raw != None:
        frame_setting = adata.raw

    return adata, frame_setting

def get_h5_genes(h5file):

    _, frame_setting = prep_h5(h5file)

    #Gets a list of all genes within the H5 files
    all_genes = frame_setting.var_names.to_list()

    return all_genes

def print_genes(genelist):

    with open(f'ALL_H5_GENES.txt', 'w+') as outfile:
        for i in genelist:
            outfile.write(f'{i}\n')
            
if __name__ == '__main__':
    args = get_args()
    allgenes = get_h5_genes(args.h5file)
    print_genes(allgenes)