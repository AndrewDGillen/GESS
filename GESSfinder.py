import csv
import math
import os
 
from collections import defaultdict

import GESS_core as GESS

from Utils.GESSargs import get_gess_args

#Gets all defined terms to skip
def parse_unacceptable(file_loc):
    dontuse = []

    with open(file_loc) as unacceptablefile:
        for line in unacceptablefile:
            term = line.strip('\n').strip(',')
            
            dontuse.append(term)
    
    return dontuse

#Runs some basic checks for arguments being in the correct format before launching into GESS in earnest
def gess_setup_check(targetannos, querydata, targetdata, h5mode, annotations):

    if targetannos == []:
        print('--GESS ERROR--\nAnnotations must be directly defined based on your supplied data')
        exit()
    
    if '.h5' in querydata or '.h5' in targetdata:

        if h5mode == None:
            print('--GESS ERROR--\nPlease define a MODE for handling data from the H5 file. GESS can report based on either:\n\tprevalence: The number of individual cells expressing a gene\n\texpression:The average expression of a gene in a given cell type')
            exit()
        
    else:
        if annotations == '':
            print('--GESS ERROR--\nAnnotations must be supplied in a seperate .csv if running from a flat expression matrix file')
            exit()

#This function allows the user to get a GESS for a given set of genes using their data using DEFAULT SETTINGS
#GESS is simply printed to terminal
def find_gess(query_gene, query_data, target_gene, target_data, annos='', unacceptable='', targetannos=[], analysis_mode=None, umithresh=1,  q_species='', t_species='', verbosity=False ):

    gess_setup_check(targetannos, query_data, target_data, analysis_mode, annos)

    if unacceptable == '':
        unacceptable_terms = ['N/A', 'NA', '', 'Mixed', 'MINIMUM', 'Many','Whole body', 'Whole']
    else:
        unacceptable_terms = parse_unacceptable(unacceptable)

    TEST = GESS.Comparison(query_gene, query_data, target_gene, target_data, unacceptable_terms, analysis_mode, umithresh, annotation_file=annos, query_species = q_species, target_species=t_species, anno_levels=targetannos)

    if not verbosity:
        gess= TEST.standard_gess(weighting=True, loglevel=2, posdistance=True, verbose=False)
       
    else:
        gess = TEST.standard_gess(weighting=True, loglevel=2, posdistance=True, verbose=True)
    
    return gess

#This function allows the user to run a test for GESS on their data, returning the GESS scores under various parameterisations to assist in choosing the correct settings.
def test_gess(query_gene, query_data, target_gene, target_data, annos, q_species, t_species, unacceptable, targetannos, analysis_mode, umithresh):

    testfilename = f'FGStests/{query_gene}_vs_{target_gene}_testdata.csv'
    if os.path.isfile(testfilename):
        os.remove(testfilename)

    if unacceptable == '':
        unacceptable_terms = ['N/A', 'NA', '', 'Mixed', 'MINIMUM', 'Many','Whole body', 'Whole']
    else:
        unacceptable_terms = parse_unacceptable(unacceptable)

    TEST = GESS.Comparison(query_gene, query_data, target_gene, target_data,unacceptable_terms, analysis_mode, umithresh, annotation_file=annos, query_species = q_species, target_species=t_species, testfile=testfilename, anno_levels=targetannos)

    for weighting in [True, False]:
        for loglevel in [2, 10]:
                for use_posdist in [True, False]:
                    TEST.test_gess(weighting, loglevel, use_posdist)

if __name__ == '__main__':
    args = get_gess_args()
            
    if args.testing:
        test_gess(args.querygene, args.querydata, args.targetgene, args.targetdata,args.annotations, args.queryspecies, args.targetspecies, args.unaceptable, args.targetannos, args.h5mode, args.umithreshold)
    else:
        gess = find_gess(
            args.querygene, 
            args.querydata, 
            args.targetgene, 
            args.targetdata, 
            args.annotations, 
            args.unacceptable, 
            args.targetannos,
            args.h5mode,
            args.umithreshold, 
            args.queryspecies, 
            args.targetspecies, 
            args.verbose)
        
        print('GESS SCORE:', f'{round(gess, 2)}%')