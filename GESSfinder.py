import csv
import math
import os
 
from collections import defaultdict
import argparse
import GESS_core_Mar24UPDATE as GESS

#Takes user arguments. 
def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-qg', '--querygene', type = str, required = True, help = 'Query Gene')
    parser.add_argument('-tg', '--targetgene', type = str, required = True, help = 'Target Gene')
    parser.add_argument('-qd', '--querydata', type=str, required=True, help ='Enrichment data corresponding to the Query genes')
    parser.add_argument('-td', '--targetdata', type=str, required=True, help ='Enrichment data corresponding to the Target genes')
    parser.add_argument('-qs', '--queryspecies', type=str, default = '', help ='The species from which query genes are derived. This only needs to be defined if not a default FlyGene species')
    parser.add_argument('-ts', '--targetspecies', type=str, default = '', help ='The species from which target genes are derived. This only needs to be defined if not a default FlyGene species')
    parser.add_argument('-a', '--annotations', type=str, default='', help = 'Tissue annotation file. Required if using a flat expression matrix')
    parser.add_argument('-test', '--testing', type=bool, default=False, help = 'Boolean switch to activate GESS testing suite.')
    parser.add_argument('-v', '--verbose', type = bool, default = False, help = 'Boolean switch activating VERBOSE mode')
    parser.add_argument('-u', '--unacceptable', type=str, default='', help="file containing uninterpretable annotations. Line seperated text file (see example)")
    parser.add_argument('-anno', '--targetannos', nargs='*', default = [], help = "Number of annotation levels")
    parser.add_argument('-h5mode', '--h5mode', type=str, choices = [None, 'prevalence','expression'], default = None, help = 'Mode of h5ad analysis - either "prevalence" (the % of reads from a context with detected expression UMI>1) or "expression" (the average expression within a given context)')
    parser.add_argument('-umi', '--umithreshold', type = int, default = 1, help = 'Sensitivity threshold for H5 parsing in expression mode. Any cell expressing a gene >= this threshold will count as an "Expressing" cell')
    args = parser.parse_args()

    return args

#Gets all defined terms to skip
def parse_unacceptable(file_loc):
    dontuse = []

    with open(file_loc) as unacceptablefile:
        for line in unacceptablefile:
            term = line.strip('\n').strip(',')
            
            dontuse.append(term)
    
    return dontuse

#This function allows the user to get a GESS for a given set of genes using their data using DEFAULT SETTINGS
#GESS is simply printed to terminal
def run_gess(q_gene, q_data, t_gene, t_data, annos, q_species, t_species, verbosity, unacceptable, levels, h5mode, umithresh):

    if unacceptable == '':
        unacceptable_terms = ['N/A', 'NA', '', 'Mixed', 'MINIMUM', 'Many','Whole body', 'Whole']
    else:
        unacceptable_terms = parse_unacceptable(unacceptable)

    TEST = GESS.Comparison(q_gene, q_data, t_gene, t_data, unacceptable_terms, h5mode, umithresh, annotation_file=annos, query_species = q_species, target_species=t_species, anno_levels=levels)

    print()
    if not verbosity:
        print('GESS SCORE:', f'{round(TEST.standard_gess(weighting=True, loglevel=2, posdistance=True, verbose=False), 2)}%')
    else:
        gess = TEST.standard_gess(weighting=True, loglevel=2, posdistance=True, verbose=True)
        print('GESS SCORE:', f'{round(gess, 2)}%')

#This function allows the user to run a test for GESS on their data, returning the GESS scores under various parameterisations to assist in choosing the correct settings.
def test_gess(q_gene, q_data, t_gene, t_data, annos, q_species, t_species, unacceptable, levels, h5mode, umithresh):

    testfilename = f'FGStests/{q_gene}_vs_{t_gene}_testdata.csv'
    if os.path.isfile(testfilename):
        os.remove(testfilename)

    if unacceptable == '':
        unacceptable_terms = ['N/A', 'NA', '', 'Mixed', 'MINIMUM', 'Many','Whole body', 'Whole']
    else:
        unacceptable_terms = parse_unacceptable(unacceptable)

    TEST = GESS.Comparison(q_gene, q_data, t_gene, t_data,unacceptable_terms, h5mode, umithresh, annotation_file=annos, query_species = q_species, target_species=t_species, testfile=testfilename, anno_levels=levels)

    for weighting in [True, False]:
        for loglevel in [2, 10]:
                for use_posdist in [True, False]:
                    TEST.test_gess(weighting, loglevel, use_posdist)

if __name__ == '__main__':
    args = get_args()

    if args.targetannos == []:
        print('--GESS ERROR--\nAnnotations must be directly defined based on your supplied data')
        exit()
    
    if '.h5' in args.querydata or '.h5' in args.targetdata:

        if args.h5mode == None:
            print('--GESS ERROR--\nPlease define a MODE for handling data from the H5 file. GESS can report based on either:\n\tprevalence: The number of individual cells expressing a gene\n\texpression:The average expression of a gene in a given cell type')
            exit()
        
    else:
        if args.annotations == '':
            print('--GESS ERROR--\nAnnotations must be supplied in a seperate .csv if running from a flat expression matrix file')
            exit()
            
    if args.testing:
        test_gess(args.querygene, args.querydata, args.targetgene, args.targetdata,args.annotations, args.queryspecies, args.targetspecies, args.unaceptable, args.targetannos, args.h5mode, args.umithreshold)
    else:
        run_gess(args.querygene, args.querydata, args.targetgene, args.targetdata,args.annotations, args.queryspecies, args.targetspecies, args.verbose, args.unacceptable, args.targetannos, args.h5mode, args.umithreshold)
