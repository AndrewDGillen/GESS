import numpy as np
import matplotlib.pyplot as plt
import scipy
import argparse

#This is optional but speeds up clustering for large matrices
import fastcluster

import pandas as pd
import seaborn as sns
import gc

from tqdm import tqdm

import GESS_core as GESS
import Utils.HandleH5s as GESS_h5

from Utils.GESSargs import get_gess_args

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
        query_gene = i_gene
        query_label = i_gene
        query_labels.append(query_label)

        for i_target in target_genes:
            
            #Target (gene, label) tuples are parsed.
            #If the target labels have not been parsed, we add these to the target label list
            target_gene = i_target
            target_label = i_target

            if tlabels_popped == False:
                target_labels.append(target_label)
            
            #Uses the GESS core to carry out individual comparisons between query and target genes
            #A Comparison object is instantiated, then the get_flygene score method is called
            #try:
            compare_object = GESS.Comparison(query_gene, args.querydata, target_gene, targetdata, args.unacceptable, args.analysis_mode, args.umithreshold, args.annotations, args.queryspecies, args.targetspecies, anno_levels=args.targetannos)
            
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
def matricise_all(query_genes, target_genes, target_data, args):
    
    #The user gene information is parsed, valid genes are identified, and only valid genes are used to generate the matrix 
    #This process is carried out for both query and target genes
    print('\nChecking your Genes of Interest')

    if '.h5' in args.querydata:
        valid_genes = GESS_h5.get_h5_genes(args.querydata)
        q_valid = list(set(query_genes).intersection(set(valid_genes)))
    else:
        q_valid = GESS.get_valid_genes(args.querydata)
        query_genes = [x for x in query_genes if x in q_valid]

    if query_genes == []:
        print('\n~~~~~ERROR~~~~~\n')
        print('No valid QUERY genes found. Please check the format of your input list, and that these genes can be found in your QUERY data file')
        #exit()

    if '.h5' in target_data:
        valid_genes = GESS_h5.get_h5_genes(target_data)
        q_valid = list(set(target_genes).intersection(set(valid_genes)))
    else:
        t_valid = GESS.get_valid_genes(target_data)
        target_genes = [x for x in target_genes if x in t_valid]

    if target_genes == []:
        print('\n~~~~~ERROR~~~~~\n')
        print('No valid TARGET genes found. Please check the format of your input list, and that these genes can be found in your TARGET data file')
        #exit()

    total_valid = len(set(query_genes).union(set(target_genes)))

    #Creates a clustered heatmp of GESS scores in query vs target lists using the Seaborn CLUSTERMAP function
    get_matrix(query_genes, target_genes, target_data, args)

    #If a filename has been supplied, we use this to save the generated matrix
    if args.save != '':
        plt.savefig(args.save)
    else:
        #Displays the generated matrix
        plt.show()

#A function to quickly check that the user settings appear to be valid for matrix generation
def test_validity(args):

    #Checks that the user has included "Query" data
    if args.input_query == '' or args.querydata == '':
        print('\n~~~~~ERROR~~~~~\n')
        print('You must supply both query genes (as a line seperated .txt file or as a python list) AND query data(as .csv or .h5ad)')
        exit()

    #This block of code deals specifically with identifying whether one or two data sets are needed.
    #If no seperate target list is provided, the matrix will be formed all-against-all from query genes
    if args.input_target == '':
        target_in = args.input_query
        target_data = args.querydata
        
    #If a seperate target list IS provided, the matrix will be all queries vs all targets
    else:
        target_in = args.input_target

        if args.targetdata != '':
            target_data = args.targetdata
        else:
            target_data = args.querydata
        
    return target_in, target_data

#The function for calling matricise from python. 
#Just builds an argparse object from user arguments and passes it on to matricise_all
def gess_matricise(query_genes,querydata, target_genes='', target_data='', queryspecies='', targetspecies='', annotations='',unacceptable='', targetannos=[], analysis_mode=None,umi_threshold=1, savefilename='', cluster_method='weighted',labelling=False):
    args = argparse.ArgumentParser()

    if type(query_genes) == str:
        query_genes = get_list_from_file(query_genes) 
    
    args.input_query = query_genes
    args.querydata = querydata

    args.input_target = target_genes
    args.targetdata = target_data

    target_in, target_data = test_validity(args)

    if type(target_in) == str:
        target_genes = get_list_from_file(target_in)
    else:
        target_genes = target_in

    args.queryspecies = queryspecies
    args.targetspecies = targetspecies
    args.annotations = annotations
    args.unacceptable = unacceptable
    args.targetannos = targetannos
    args.analysis_mode = analysis_mode
    args.umithreshold = umi_threshold
    args.save = savefilename
    args.cluster_method = cluster_method
    args.label = labelling

    matricise_all(query_genes, target_genes, target_data, args)

if __name__=='__main__':

    args=get_gess_args()

    target_in, target_data = test_validity()

    query_genes = get_list_from_file(args.input_query)
    target_genes = get_list_from_file(args.target_in)

    matricise_all(query_genes, target_genes, target_data, args)
