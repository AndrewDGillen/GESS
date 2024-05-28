import numpy as np
import matplotlib.pyplot as plt
import scipy
import argparse

import GESS_core_Mar24UPDATE as GESS
import GESS_h5handling_test as GESS_h5
import pandas as pd
import seaborn as sns
import gc

from tqdm import tqdm

#This is optional but speeds up clustering for large matrices
import fastcluster

#Takes user arguments
def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-iq', '--input_query', type=str, required=True, help = 'Chosen genes from the query species you wish to include within the matrix')
    parser.add_argument('-it', '--input_target', type=str, default='', help = 'Chosen genes from the target species you wish to include within the matrix')
    parser.add_argument('-qd', '--querydata', type=str, required=True, help ='Enrichment data corresponding to the Query genes')
    parser.add_argument('-td', '--targetdata', type=str, default = '', help ='Enrichment data corresponding to the Target genes')
    parser.add_argument('-a', '--annotations', type=str, default='FlyGene_Annotations.csv', help = 'Tissue annotation file. Uses FlyGene tissue annotations by default')
    parser.add_argument('-qs', '--queryspecies', type=str, default = '', help ='The species from which query genes are derived. This only needs to be defined if not a default FlyGene species')
    parser.add_argument('-ts', '--targetspecies', type=str, default = '', help ='The species from which target genes are derived. This only needs to be defined if not a default FlyGene species')
    parser.add_argument('-s', '--save', type = str, default='', help = 'Filename under which to save the produced matrix. Please only use image file formats. Leave blank if you do not wish to save')
    parser.add_argument('-c', '--cluster_method', type = str, choices = ['single', 'complete', 'median', 'average', 'weighted', 'centroid', 'ward'], default = 'weighted', help = 'Clustering method for heatmap/dendogram generation')
    parser.add_argument('-l', '--label', type=bool, default=False, help = 'Boolean switch for activating labelling of the heatmap. This gets very busy for big heatmaps!')
    parser.add_argument('-u', '--unacceptable', type=str, default='', help="file containing uninterpretable annotations. Line seperated text file (see example)")
    parser.add_argument('-n', '--numberannos', type = int, default = 3, help = "Number of annotation levels")
    parser.add_argument('-anno', '--targetannotations', nargs='*', default = [], help = "h5 attributes to use as annotation levels for the data. This is REQUIRED if using a .h5 file as input")
    parser.add_argument('-h5mode', '--h5mode', type=str, choices = [None, 'prevalence','expression'], default = None, help = 'Mode of h5ad analysis - either "prevalence" (the % of reads from a context with detected expression UMI>1) or "expression" (the average expression within a given context)')
    parser.add_argument('-umi', '--umithreshold', type = int, default = 1, help = 'Sensitivity threshold for H5 parsing in expression mode. Any cell expressing a gene >= this threshold will count as an "Expressing" cell')
 
    return parser.parse_args()

#Given user input data, we calculate the GESS between each query gene and each target gene.
#GESS scores are used to calculate euclidean distance in expression pattern similarity between scores
#Genes are then clustered based on their GESS to every other gene, allowing a clustermap to be printed which represents the similarity between expression patters across the two sets of genes.
def get_matrix(query_genes, target_genes, targetdata, args):
    #A list is prepared which will hold all GESSs for all query genes
    all_comparisons = []

    #We initialise lists of labels, and create a boolean flag for ppopulation of target labels
    query_labels = []
    target_labels = []
    tlabels_popped = False

    #For each query gene, we calculate GESS to each target gene
    #Note that I use TQDM here for a progress bar - this is purely preference and could be removed without issue
    for i_gene in tqdm(query_genes):

        #A list is set up to hold all the GESS for a given gene
        all_gene_comparisons = []

        #The (gene, label) tuple is parsed and label is stored
        query_gene = i_gene[0]
        query_label = i_gene[1]
        query_labels.append(query_label)

        for i_target in target_genes:
            
            #Target (gene, label) tuples are parsed.
            #If the target labels have not been parsed, we add these to the target label list
            target_gene = i_target[0]
            target_label = i_target[1]

            if tlabels_popped == False:
                target_labels.append(target_label)
            
            #Uses the GESS core to carry out individual comparisons between query and target genes
            #A Comparison object is instantiated, then the get_flygene score method is called
            #try:
            compare_object = GESS.Comparison(query_gene, args.querydata, target_gene, targetdata, args.unacceptable, args.h5mode, args.umithreshold, args.annotations, args.queryspecies, args.targetspecies, anno_levels=args.targetannotations)
            
            try:
                gess = compare_object.standard_gess(weighting=True, loglevel=2, posdistance=True, verbose=False)
            except ZeroDivisionError:
                print(f'No categories in common between {i_gene} and {target_gene}')
                gess = 0

            
            gc.collect()

            all_gene_comparisons.append(gess)

        #After a single query gene, all target labels will be gathered. As such, we no longer add to the target label list
        tlabels_popped = True

        all_comparisons.append(all_gene_comparisons)

    #The all comparisons list is converted to a numpy array
    harvest = np.array(all_comparisons)

    print(harvest)
    #Generates individual linkage matrices to allow optimal ordering (based on euclidean distance) for the clustermap
    ordered_columns = scipy.cluster.hierarchy.linkage(np.transpose(harvest), method=args.cluster_method, metric='euclidean', optimal_ordering=True)
    ordered_rows = scipy.cluster.hierarchy.linkage(harvest,method=args.cluster_method, metric='euclidean', optimal_ordering=True)
    
    #Generates a pandas datafame containing the gess scores, labelled columns and rows as appropriate
    df = pd.DataFrame(harvest, columns = target_labels)
    df.index=query_labels

    #Generates a Seaborn Clustermap from the supplied data. Parameters:
        #Colour range 0-100 (as GESS are percentages)
        #Labels on each row and column
        #Clustering for both rows and columns based on the user-selected clustering algorith, 
        #   by default we use the Weighted clustering algorithm/ WPGMA (Weighted Pair Group Method with Arithmetic Mean) 
        #Uses pre-generated linkage matrices to preserve optimal (ie minimal) euclidean distance between neighbors
    
    #Note that this clustering CAN use fastcluster, but doesn't NEED to
    sns.clustermap(df, vmin=0, vmax=100, annot=args.label, method=args.cluster_method, xticklabels = 1, yticklabels=1, row_cluster=True, col_cluster=True, row_linkage = ordered_rows, col_linkage=ordered_columns)

#Parses a 'Matricise-format' .txt file into usable (gene, label) tuples.
#These Files MUST contain one gene/line, and MAY OPTIONALLY contain A SINGLE "Label" for a gene (which can be any identifier) following A SINGLE \t
#If no labels are provided, genes are just labeled by their ID
def get_list_from_file(targetfile):
    resulting_list = []

    with open(targetfile) as filein:
        for line in filein:
            cleanline = line.strip('\n')

            if '\t' in cleanline:
                geneandlabel = (cleanline.split('\t')[0],cleanline.split('\t')[1])
            else:
                geneandlabel = (cleanline, cleanline)

            resulting_list.append(geneandlabel)
    
    return resulting_list

#Handles generating and saving a matrix from the user supplied data
def matricise_all(target_in, target_data, args):
    
    #The user gene information is parsed, valid genes are identified, and only valid genes are used to generate the matrix 
    #This process is carried out for both query and target genes
    print('\nChecking your Genes of Interest')

    query_genes = get_list_from_file(args.input_query)

    if '.h5' in args.querydata:
        valid_genes = GESS_h5.get_h5_genes(args.querydata)
        q_valid = list(set(query_genes).intersection(set(valid_genes)))
    else:
        q_valid = GESS.get_valid_genes(args.querydata)
        query_genes = [x for x in query_genes if x[0] in q_valid]

    if query_genes == []:
        print('\n~~~~~ERROR~~~~~\n')
        print('No valid QUERY genes found. Please check the format of your input list, and that these genes can be found in your QUERY data file')
        exit()

    target_genes = get_list_from_file(target_in)

    if '.h5' in target_data:
        valid_genes = GESS_h5.get_h5_genes(target_data)
        q_valid = list(set(target_genes).intersection(set(valid_genes)))
    else:
        t_valid = GESS.get_valid_genes(target_data)
        target_genes = [x for x in target_genes if x[0] in t_valid]

    if target_genes == []:
        print('\n~~~~~ERROR~~~~~\n')
        print('No valid TARGET genes found. Please check the format of your input list, and that these genes can be found in your TARGET data file')
        exit()

    total_valid = len(set(query_genes).union(set(target_genes)))
    print(f'Genes of Interest Checked. {total_valid} are valid in these datasets\n')

    #Creates a clustered heatmp of GESS scores in query vs target lists using the Seaborn CLUSTERMAP function
    get_matrix(query_genes, target_genes, target_data, args)

    #If a filename has been supplied, we use this to save the generated matrix
    if args.save != '':
        plt.savefig(args.save)
    else:
        #Displays the generated matrix
        plt.show()

if __name__=='__main__':

    args=get_args()

    #Warns the user about how long H5 processing can take
    if '.h5' in args.querydata or '.h5' in args.targetdata:
        print("\n--WARNING--\nProcessing large .h5 files can be very time consuming, especially if you have many genes of interest.\nThis can take up to 1 minute PER COMPARISON \nIf this becomes an issue, consider processing these in smaller chunks.")
    
    #This block of code deals specifically with identifying whether one or two data sets are needed.
    #If no seperate target list is provided, the matrix will be formed all-against-all from query genes
    if args.input_target == '':
        target_in = args.input_query
        target_data = args.querydata
    
    #If a seperate target list IS provided, the matrix will be all queries vs all targets
    else:
        target_in = args.input_target
        target_data = args.targetdata
        
        #Naturally, this requires target data be supplied. If none is, the program will exit out
        if target_data == '':
            print('\n~~~~~ERROR~~~~~\n')
            print('If you wish to use seperate target genes, please supply a target data file using -t')
            exit()

    matricise_all(target_in, target_data, args)
