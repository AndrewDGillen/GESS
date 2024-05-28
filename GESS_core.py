import csv
import math
import os
import sys

from collections import defaultdict
from utils.GESS_HandleH5s.py import H5_GESS
import pandas as pd

#A class for handling GESS calculations
class Comparison():

    def __init__(self, query_gene, query_data,  target_gene, target_data, unnaceptable_terms, h5mode = None, umi_thresh = 1, annotation_file='', query_species = '', target_species='',testfile = '', anno_levels = ['Name']):
        
        if annotation_file != '':
            self.annotation_file = annotation_file
        else:
            self.annotation_file ='FlyGene_Annotations.csv'

        self.unnaceptable_terms = unnaceptable_terms
        self.anno_levels=anno_levels

        self.testfile = testfile
        self.querydata = query_data
        self.targetdata = target_data
        self.h5mode = h5mode

        #Genes to compare - a "Query" gene and a "Target" gene. Effectively these are equivalent, and distinction is just for clarity
        self.querygene = query_gene
        self.targetgene = target_gene

        #If the user fails to supply them, the gene name format is used to infer species
        #In cases where this is not baked into FlyGene, the program will error out.

        if query_species == '' and self.h5mode == None:
            self.queryspecies= self.get_species(self.querygene)
        else:
            self.queryspecies = query_species
        
        if target_species == '' and self.h5mode == None:
            self.targetspecies = self.get_species(self.targetgene)
        else:
            self.targetspecies = target_species

        #Handles datafiles differently for h5ad files vs flat expression matrices
        if '.h5' in self.querydata:
            query_alllevels = H5_GESS(anno_levels, self.querydata, self.querygene, self.h5mode, umi_thresh)

        else:
            #Populates alternative annotations for tissues in a query species
            #By default, these are pulled from the FlyGene annotations, but a user can also supply their own.
            self.queryannotations, _ = populate_annotations(self.queryspecies, self.targetspecies, self.annotation_file, self.unnaceptable_terms, self.anno_levels)

            if len(self.queryannotations) < 1:
                print('--GESS ERROR--\nNo annotations could be identified for the query dataset from the supplied annotation file\nPlease be very careful with column name/formatting, and ensure entries are in the same format as the supplies argument')
                exit()

            #vector of enrichments from individual tissues in format (Category name, enrichment)
            #From QUERY species tissues
            #Query minimum is currently set to 0 in all cases
            self.querytissues = file_to_vectors(self.querydata,self.queryannotations, self.querygene)

            #Data is sorted into other vectors
            query_alllevels = tissues_to_annotations(self.querytissues, self.queryannotations, self.anno_levels, self.unnaceptable_terms)

        #Similar variables are assembled for the target sample
        if '.h5' in self.targetdata:
            target_alllevels = H5_GESS(anno_levels, self.targetdata, self.targetgene, h5mode, umi_thresh)

        else:
            #Populates alternative annotations for tissues in target species
            #By default, these are pulled from the FlyGene annotations, but a user can also supply their own.
            _, self.targetannotations = populate_annotations(self.queryspecies, self.targetspecies, self.annotation_file, self.unnaceptable_terms, self.anno_levels)

            #From TARGET species tissues
            #Query minimum is currently set to 0 in all cases
            self.targettissues = file_to_vectors(self.targetdata, self.targetannotations, self.targetgene)
            
            #Data is sorted into other vectors - one for Function, one for Type
            target_alllevels = tissues_to_annotations(self.targettissues, self.targetannotations,  self.anno_levels, self.unnaceptable_terms)

        self.all_comparisons = []
        #Each comparable pair (Tissue, Function and Type) is assembled into a tuple of (query, target)
        #These tuples are then collected into a list

        for i in anno_levels :
            self.all_comparisons.append((query_alllevels[i], target_alllevels[i]))#


    #FOR TESTING - 
    # #Allows manual setting of data vectors based on data from files
    # def update_datavectors(self, querymatrix, targetmatrix): 
    #     self.querytissues, self.queryfuncs, self.querytypes = querymatrix
    #     self.targettissues, self.targetfuncs, self.targettypes = targetmatrix

    #     self.tissues_compare = (self.querytissues, self.targettissues)
    #     self.func_compare = (self.queryfuncs, self.targetfuncs)
    #     self.type_compare = (self.querytypes, self.targettypes)

    #     self.all_comparisons = [(self.tissues_compare, False), (self.func_compare, True), (self.type_compare,True)]
    
    #Uses gene name to infer species in known cases. 
    #At present, we only have D. mel and A. aeg baked in - these are the species currently worked up for FlyGene
    def get_species(self, genename):
        geneformats = {
            #'FBGN':('Drosophila melanogaster', 'Enrichment_FPKM_Gene_Female.csv'),
            'FBGN':('Drosophila melanogaster'),
            'AAEL':('Aedes aegypti')
            }

        try:
            species_name = geneformats[genename[:4].upper()]

        except:
            print('\n~~~ERROR~~~')
            print('The gene format cannot be assigned to a FlyGene species')
            print('Please manually provide species names if not using a FlyGene species')

        return species_name
    
    #********
    #TBD - Update definitiion
    #********

    #FlyGene Score is defined as:
    #       1/ 1 + (|Position of category in Query gene data - Position of category in Target gene data|)
    #       PLUS 1/ 1 + ABSOLUTE DIFFERENCE(log2(Enrichment of category in Query gene data) - log2(Enrichment of category in Target gene data))
    #       FOR EACH CATEGORY (That is, EVERY tissue, EVERY Functional annotation, EVERY tissue type) present in both query and target category sets
    #       FOR EACH ANNOTATION LEVEL (That is, all of Tissue, Function and Type)
    
    #Note that we currently log2 transform enrichment values for this calculation. This achieves two things:
        #Corrects the calculation to represent the DISTANCE of both enrichment and depletion from 1 (ie no change from whole insect)
        #Fixes the scale to fully represent depletion instead of it being bounded between 0-1

    #This produces a score which is proportional to the similarity in gene expression patterns.
    #Fundamentally, the upper limits of this score are limited by the number of shared categories, and therefore by the depth of each dataset.
    #The upper limit for fly gene score is defined:
    #       (2 * shared categories at a given annotation level)
    #       SUMMED ACROSS ANNOTATION LEVELS

    #FlyGene Score is currently presented as a percentage of this maximum BUT this may overrepresent similarity in expression pattern between shallower data sets.

    #A special version of GESS is provided for testing, which A - gives very verbose printouts to terminal and B - Prints test results to file
    def test_gess(self, weighting, loglevel, posdistance):
        print(f'{self.querygene} vs {self.targetgene}')
        
        with open(self.testfile, 'a') as test_file:
            INtest = csv.writer(test_file)

            INtest.writerow(['TEST PARAMETERS', f'Weighting:   {weighting}', f'Log base {loglevel}', f'Using Position Distancing:   {posdistance}',])
            INtest.writerow([])

            #We initialise score summaries for both actual and maximum genescore 
            full_genescore_add = 0
            full_genescore_multiplied = 0

            max_genescore_add = 0
            max_genescore_multiply = 0

            #Individual annotation levels are processed seperately, allowing scores to be summed independently
            for annotation_level in self.all_comparisons:

                comparison_pair = annotation_level

                #get_common_catlist is called to find the intersection between the categories available for the query and target genes
                consistent_cats = get_common_catlist(comparison_pair)
                print('\n\n')
                print('----Comparing categories:', consistent_cats)
                
                #For annotation levels where multiple tissues may fall under the same category (IE Function or Type level comparisons), we gather AVERAGE ENRICHMENT for each category
                comparison_pair = (condense_cats(comparison_pair[0]), condense_cats(comparison_pair[1]))

                #We remove all non-matched categories from both query and target lists
                #Leaving these categories in the list would affect positioning and thus invalidate comparisons.
                scrubbed_queries, scrubbed_targets = scrub_lists(comparison_pair, consistent_cats)

                #The summed annotation level FlyGene Score component is gathered by comparing the query vs target vectors
                comparison_score_sum, comparison_score_factor, usedcats = make_comparison(scrubbed_queries, scrubbed_targets, weighting, loglevel, posdistance, self.h5mode, INtest)

                #The maximum genescore for a comparison is equal to 2* the number of categories available for comparison
                max_genescore_add += 2* (len(usedcats))

                #When multiplying the Position and Enrichment scores, the max score is simply the number of categories
                max_genescore_multiply += len(usedcats)

                #The full genescore is incremented by this summed comparison score
                full_genescore_add += comparison_score_sum
                full_genescore_multiplied += comparison_score_factor
                INtest.writerow([])
            
            INtest.writerow(['FULL ADDITION SCORE:', full_genescore_add, 'MAX ADDITION SCORE:', max_genescore_add, 'PERCENTAGE SCORE:', (full_genescore_add/max_genescore_add) * 100])
            print('\n\n\n')
            print('****ADDING*****')
            print('FGS:', full_genescore_add)
            print('THEORETICAL MAXIMUM: \t', max_genescore_add)
            print('PERCENTAGE: ', (full_genescore_add/max_genescore_add) * 100)

            INtest.writerow(['FULL MULTIPLICATIVE SCORE:', full_genescore_multiplied, 'MAX ADDITION SCORE:', max_genescore_multiply, 'PERCENTAGE SCORE:', (full_genescore_multiplied/max_genescore_multiply) * 100])
            print('\n****MULTIPLICATION*****')
            print('FGS:', full_genescore_multiplied)
            print('THEORETICAL MAXIMUM : \t', max_genescore_multiply)
            print('PERCENTAGE: ', (full_genescore_multiplied/max_genescore_multiply) * 100)

            for i in range(3):
                INtest.writerow([])

            return (full_genescore_multiplied/max_genescore_multiply) * 100

    #GESS for non-test cases - no printing to various files, only returning the calculated GESS using the provided settings
    def standard_gess(self, weighting, loglevel, posdistance, verbose):

        #We initialise score summaries for both actual and maximum genescore 
        full_genescore_multiplied = 0
        max_genescore_multiply = 0

        #Individual annotation levels are processed seperately, allowing scores to be summed independently
        for annotation_level in self.all_comparisons:

            comparison_pair = annotation_level

            #get_common_catlist is called to find the intersection between the categories available for the query and target genes
            consistent_cats = get_common_catlist(comparison_pair)
            
            #For annotation levels where multiple tissues may fall under the same category (IE Function or Type level comparisons), we gather AVERAGE ENRICHMENT for each category
            comparison_pair = (condense_cats(comparison_pair[0]), condense_cats(comparison_pair[1]))
            # print(annotation_level)
            # print()
            # print(comparison_pair) 
            # print()
            # print(consistent_cats)
            #We remove all non-matched categories from both query and target lists
            #Leaving these categories in the list would affect positioning and thus invalidate comparisons.
            scrubbed_queries, scrubbed_targets = scrub_lists(comparison_pair, consistent_cats)
            
            if verbose == False:
                blockPrint()
                #pass
            else:
                print('\n\n')
                print('----Comparing categories:', consistent_cats)

            #The summed annotation level FlyGene Score component is gathered by comparing the query vs target vectors
            comparison_score_sum, comparison_score_factor, usedcats = make_comparison(scrubbed_queries, scrubbed_targets, weighting, loglevel, posdistance, self.h5mode)

            #When multiplying the Position and Enrichment scores, the max score is simply the number of categories used
            max_genescore_multiply += len(usedcats)
            
            enablePrint()

            #The full genescore is incremented by this summed comparison score
            full_genescore_multiplied += comparison_score_factor

        return (full_genescore_multiplied/max_genescore_multiply) * 100

#Populates enrichment data from supplied data files by annotation
def file_to_vectors(testdata, annotations_by_tissue, geneofinterest):
    tissue_by_index = {}
    tissuelist = []

    with open(testdata) as testfile:
        
        if '.csv' in testdata:
            INdata = csv.reader(testfile)
        elif '.tsv' in testdata:
            INdata=csv.reader(testfile, delimiter='\t')
        else:
            print('---ERROR---')
            print('Data filetype not supported')
            exit()

        linecount = 0
                
        for line in INdata:

            if linecount == 0:

                itemcount = 0

                for item in line: 
                    tissue_by_index[itemcount] = item
                    
                    itemcount += 1

                    
                linecount += 1


            elif line[0] == geneofinterest:

                itemcount = 1

                for item in line[1:]:
                    
                    try:
                        tissue = tissue_by_index[itemcount]

                        enrichment = float(item)
                            
                        tissuelist.append((tissue, enrichment))
                        
                        itemcount += 1
                    except:
                        pass

    return sort_generic_list(tissuelist)

# Disable printing
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore printing
def enablePrint():
    sys.stdout = sys.__stdout__

#A very simple function to sort a provided list of (Category, Enrichment) tuples based on enrichment values, with highest values brought to the top
def sort_generic_list(incoming_list):
    return sorted(incoming_list, key=lambda x:x[1], reverse=True)
        
#Finds the intersection between two sets of categories - ie, the tissues/functions/types shared between query and target sets
def get_common_catlist(listpair):
    querylist, targetlist = listpair

    all_query_cats = set(x[0] for x in querylist)
    all_target_cats = set(x[0] for x in targetlist)

    return all_query_cats.intersection(all_target_cats)

#given a list of (Category, Enrichment) tuples in which several categories may be identical, this function averages the enrichment for each category and returns a sorted list of tuples
def condense_cats(individual_list):

    #We create a dictionary of summed enrichment values by category (condensation)
    #As well as a dictionary containing the number of entries by category (number_of_entries
    condensation = defaultdict(lambda:0)
    number_of_entries = defaultdict(lambda:0)

    #These dictionaries are populated using the provided list
    for category in individual_list:
        condensation[category[0]] += category[1]
        number_of_entries[category[0]] += 1
    
    #We create a new list of (Category, Enrichment) tuples, in this case using AVERAGE ENRICHMENT across all entries from the initial list with matching categories
    unsorted_list = [(x, condensation[x]/number_of_entries[x]) for x in condensation]

    #We sort this list by average enrichment value and return it
    return sort_generic_list(unsorted_list)

#Given a pair of lists (query, target) and a list of matched categories, we remove all non-matched categories from each list and return the resulting list
def scrub_lists(listpair, cats_in_use):
    querylist, targetlist = listpair

    fixed_querylist = [x for x in querylist if x[0] in cats_in_use]
    fixed_targetlist = [x for x in targetlist if x[0] in cats_in_use]

    #These should not require re-sorting as the original list order is preserved, just minus all non-matching tissues
    return fixed_querylist, fixed_targetlist

def recursive_checkback(vector, position_checked):

    current_position = position_checked
    current_enrichment = vector[current_position]

    next_position = position_checked +1

    if next_position <= len(vector)-1:
        next_enrichment = vector[next_position]

        if next_enrichment == current_enrichment and next_position <= len(vector)-1:

            current_position = recursive_checkback(vector, next_position)
        
        elif next_enrichment == current_enrichment:
            current_position = next_position

    return current_position

#Gets the "Uniqueness" of a given position.
#This is defined as the MINIMUM of its distance from either the position before it or the position after it
#Where distance = higher enrichment/smaller enrichment
def get_uniqueness(vector, position):

    all_possible_distances = []

    if position > 0:
        distance_from_prev = abs(vector[position-1]-vector[position]) 
        all_possible_distances.append(distance_from_prev)
    
    if position < len(vector)-1:
        distance_from_next = abs(vector[position] - vector[position+1])
        all_possible_distances.append(distance_from_next)
    
    if all_possible_distances == []:
        return 1
    
    return min(all_possible_distances)

#This function carries out the main thrust of FlyGene Score calculation, generating a single summed value across all categories (ie Tissue, Function, Type) 
#This sum represents the similarity between BOTH relative inter-tissue gene expression relationships and per-tissue gene enrichment
# Again, this can be summarised as SUM[across all categories]((1/1+ abs(position difference) + (1/1+abs(enrichment difference))) 
def make_comparison(querylist, targetlist, weightoption, loglevel, use_posdist, exptype, file_obj = ''):

    #New feature - only the top 25 most "critical" categories for each gene are studied
    #This prevents information flooding from large scRNA-seq datasets, where MOST cell types will inevitably have relatively low expression
    #In turn, this led to genes which were very cell-type specific to different cell types having high GESS
    critical_categories = []

    if len(querylist) > 25 :
        critical_categories.extend([x[0] for x in querylist[:25]])
    else:
        critical_categories.extend([x[0] for x in querylist])

    if len(targetlist) > 25 :
        critical_categories.extend([x[0] for x in targetlist[:25]])
    else:
        critical_categories.extend([x[0] for x in targetlist])

    critical_categories = set(critical_categories)
    assert len(critical_categories) <=  50

    #First, lists of query and target (category, enrichment) tuples are seperated into individual vectors. 
    #Positions of individual categories will be identical in both cat and score vectors
    query_catvector = [x[0] for x in querylist if x[0] in critical_categories]
    query_scorevector = [x[1] for x in querylist if x[0] in critical_categories]
    target_catvector = [x[0] for x in targetlist if x[0] in critical_categories]
    target_scorevector = [x[1] for x in targetlist if x[0] in critical_categories]

    #Enrichment scores of 0.0 are problemation - they cannot be log transformed.
    #As such, these enrichment scores are set to the absolute minimum of all tissues in the comparison
    all_scores = [] + query_scorevector + target_scorevector
    all_scores = [x for x in all_scores if x != 0.0]

    #if none of the scores in the comparable tissues are > 0, we set the minimum possible to 0.001
    #Fundamentally, this is arbitrary - it could be anything, and still match across samples
    if all_scores == []:
        all_scores.append(0.001)

    minimum_score = min(all_scores)

    query_scorevector = [x if x != 0.0 else minimum_score for x in query_scorevector]
    target_scorevector = [x if x != 0.0 else minimum_score for x in target_scorevector]

    #Creates log base (loglevel) versions of the score vector 
    #ONLY USED for EXPRESSION or ENRICHMENT analyses
    if exptype != 'prevalence':
        qlog_scorevector = [math.log(x, loglevel) for x in query_scorevector ]
        tlog_scorevector = [math.log(x, loglevel) for x in target_scorevector]
    else:
        qlog_scorevector =  [] + query_scorevector
        tlog_scorevector = [] + target_scorevector

    print(qlog_scorevector)
    print(tlog_scorevector)

    #Initialises sums for scores, both additive and multipicative
    #Note that only multiplicative is typically used - additive scoring has proven to be less effective in ranking expression patterns.
    #Additive scores are maintained for testing and flexibility only
    crosscat_sum = 0
    crosscat_mult_sum = 0
    
    #Initialises output rows for test file
    catrow = ['CATEGORY:']
    qexpression_row = ['QUERY ENRICHMENT']
    qpositionrow = ['QUERY CATEGORY POSITION']
    qpos_uniquenessrow = ['QUERY CAT UNIQUENESS']
    qlogrow = ['QUERY LOG ENRICHMENT']
    texpression_row = ['TARGET ENRICHMENT']
    tpositionrow = ['TARGET CATEGORY POSITION']
    tpos_uniquenessrow = ['TARGET CAT UNIQUENESS']
    tlogrow = ['TARGET LOG ENRICHMENT']
    positiondiffrow = ['POSITION DIFFERENCE']
    worked_posrow = ['POSITION DIFFERENCE CORRECTED FOR UNIQUENESS']
    scorediffrow = ['SCORE DIFFERENCE']
    worked_scorerow=['WORKED SCORE DIFFERENCE']
    adddiffrow = ['CATEGORY SUM (ADDED)']
    multidiffrow = ['CATEGORY SUM (MULTIPLIED)']

    #Each individual category - ie tissue, function, or type - is processed seperately
    for category in query_catvector:
        catrow.append(category)
        print('\n')

        #We find the first position at which the category's enrichment value appears within the dataset, and use this as the category's position
        #This ensures cases in which multiple tissues have the same enrichment are not unfairly penalised in position scoring
        category_position_query  = recursive_checkback(query_scorevector,query_catvector.index(category))
        qpositionrow.append(category_position_query)

        #Unlogged scores are saved for testing reasons
        q_score_unlogged = query_scorevector[category_position_query]
        qexpression_row.append(q_score_unlogged)

        #Log base (loglevel) scores are used for actual enrichment comparisons
        #This ensures gene repression is on the same scale as enrichment -ie not bounded between 0-1
        category_score_query = qlog_scorevector[category_position_query]
        qlogrow.append(category_score_query)
        
        #The same metrics (corrected position, enrichment and log enrichment) are gathered for the target gene
        category_position_target = recursive_checkback(target_scorevector,target_catvector.index(category))
        tpositionrow.append(category_position_target)

        t_score_unlogged = target_scorevector[category_position_target]
        texpression_row.append(t_score_unlogged)  

        if exptype != 'prevalence':
            category_score_target = math.log(target_scorevector[category_position_target],loglevel)
        else:
            category_score_target = target_scorevector[category_position_target]

        tlogrow.append(category_score_target)

        #If the "reversal" weighting option is activated, we aim to increase the ORDER OF MAGNITUDE of the score difference by one. 
        #That is to say, we consider a difference in enrichment 10 times more important if one sample (query/target) is enriched and the other is repressed.
        if weightoption:

            #Given enrichments are log transformed, enrichment is any positive number, while repression is a negative.
            #As such, if 0 is between the enrichment in query and target samples, reversal weighting should be activated.
            if min(category_score_query, category_score_target) < 0 < max(category_score_query, category_score_target):
                reversal_weighting = 10
            else:
                reversal_weighting = 1
        else:
            reversal_weighting = 1
        
        #Score difference is equal to reversal weighting times the absolute difference between query and target log(enrichments)
        score_difference =  reversal_weighting * (abs(category_score_query-category_score_target))

        #Position difference is equal to the absolute difference betwee query and target positions for a given sample
        #Ie |Xq -Xt| where a given sample is the X most enriched in a whole set of samples
        position_difference = abs(category_position_query-category_position_target)

        #We then calculate the "uniqueness" of each sample - how different its enrichment is from the positions immediately preceeding and following it
        qposition_uniqueness = get_uniqueness(qlog_scorevector, category_position_query)
        tposition_uniqueness = get_uniqueness(tlog_scorevector, category_position_target)

        #These are added to a printable row
        qpos_uniquenessrow.append(qposition_uniqueness)
        tpos_uniquenessrow.append(tposition_uniqueness)
        
        #If position uniqueness is included in GESS calculation (it IS by default), the full POSITION function is:
        #1/1+ |position difference times (Maximum uniqueness across both samples/10)|

        #If not,  the POSITION function is:
        #1/1+position difference
        if use_posdist == True:
            max_uniqueness = max([qposition_uniqueness, tposition_uniqueness])/10
            position_function = 1/(1+abs(position_difference*max_uniqueness))
        else:
            position_function = 1/(1+position_difference)
        
        #If verbose mode is activated, relevant figures are printed to stdout
        print('POSITIONAL INFO')
        print(
            f'CATEGORY: {category}\t', f'QUERY POSITION: {category_position_query}\t', f'QUERY UNIQUENESS: {qposition_uniqueness}', f'TARGET POSITION: {category_position_target}\t', 
            f'TARGET UNIQUENESS: {tposition_uniqueness}',f'ABSOLUTE DIFFERENCE: {position_difference}\t', f'SCORED for summary: {position_function}'
            )
        
        #The score function is equal to:
        #10/10+|score difference|
        #Note that score difference has already been modified by the reversal weighting
        #Fundamentally, score difference is reduced to order of magnitude difference by using 10/10+... 
        if exptype == 'prevalence':
            score_difference *= 100
            
        score_function = 10/ (10+abs(score_difference))
        
        #relevant testing output rows are updated
        positiondiffrow.append(position_difference)
        worked_posrow.append(position_function)
        scorediffrow.append(score_difference)
        worked_scorerow.append(score_function)

        #If verbose mode is activated, relevant figures are printed to stdout
        print('ENRICHMENT INFO')
        print(
            f'CATEGORY: {category}\t', f'QUERY ENRICHMENT: {q_score_unlogged}\t', f'Log2 QUERY ENRICHMENT: {category_score_query}\t', f'TARGET ENRICHMENT: {t_score_unlogged}\t', 
            f'Log2 TARGET ENRICHMENT: {category_score_target}\t', f'ABSOLUTE DIFFERENCE: {score_difference}\t', f'SCORED for summary: {score_function}'
            )   
        
        #The category sum is equal to position funciton plus score function
        #THIS IS NOT USED AS STANDARD
        category_sum = position_function + score_function
        adddiffrow.append(category_sum)

        #Crosscat sum - ie the sum across all categories at a given annotation level - is updated by adding the ADDITIVE CATEGORY SUM
        #THIS IS NOT USED AS STANDARD
        crosscat_sum += category_sum
        
        #The multiplicative score is equal to position funciton times score function
        #THIS IS THE DEFAULT GESS METRIC USED
        category_mult = position_function * score_function
        multidiffrow.append(category_mult)

        #Crosscat MULT sum - ie the sum across all categories at a given annotation level - is updated by adding the CATEGORY MULTIPLICATIVE RESULT
        #THIS IS THE DEFAULT GESS METRIC USED
        crosscat_mult_sum += category_mult

    #If testing is enabled, the output file is updated.
    if file_obj != '':
        for list in [catrow, qexpression_row, qpositionrow, qpos_uniquenessrow, qlogrow, texpression_row, tpositionrow, tpos_uniquenessrow, tlogrow, positiondiffrow, worked_posrow, scorediffrow, worked_scorerow, adddiffrow, multidiffrow]:
            file_obj.writerow(list)
    print(crosscat_sum, crosscat_mult_sum, critical_categories)

    return crosscat_sum, crosscat_mult_sum, critical_categories

#Using the defined annotation file, we find the tissues which are present in the query species and target species
def populate_annotations(qspecies, tspecies, annotation_location, badterms, annolevels):
    
    query_annotations = defaultdict(dict)
    target_annotations = defaultdict(dict)

    annoframe = pd.read_csv(annotation_location)

    if qspecies not in annoframe.columns:
        print(f'--GESS ERROR--\nSpecies "{qspecies}" is not a column in your annotation file.')
        exit()
    if tspecies not in annoframe.columns:
        print(f'--GESS ERROR--\nSpecies "{tspecies}" is not a column in your annotation file.')
        exit()
    
    for annotation in annolevels:

        if annotation not in annoframe.columns:
            print(f'--GESS ERROR--\nAnnotation "{annotation}" is not a column in your annotation file.')
            exit()

        def type_caller(df_row):

            if df_row[qspecies].upper() == 'YES' and df_row[annotation] not in badterms:
               query_annotations[annotation][df_row['Name']]=df_row[annotation]
            
            if df_row[tspecies].upper() == 'YES':
               target_annotations[annotation][df_row['Name']]=df_row[annotation]
        
        annoframe.apply(type_caller, axis=1)

    return query_annotations, target_annotations

#With a list of tissues shared between species, and the annotations_by_tissue dictionary, we populate individual lists for shared functions and shared types.
def tissues_to_annotations(list_by_tissue, annotations_by_tissue, anno_levels, unacceptable):

    all_lists = defaultdict(list)

    all_lists[0] = [x for x in list_by_tissue if x[0] not in unacceptable]

    for entry in list_by_tissue:
        base_anno = entry[0]
        enrichment = entry[1]

        for i in anno_levels:

            try:
                appropriate_anno = annotations_by_tissue[i][base_anno]
                all_lists[i].append((appropriate_anno, enrichment))
            except:
                pass
        
    for i_list in all_lists:
        all_lists[i_list] = sort_generic_list(all_lists[i_list])

    return all_lists


#A broadly usable funciton which takes an enrichment data file, harvests all gene names, and returns a set of genes for which data exists.
#This is useful to prevent innapropriate search execution when a gene has no data associated
def get_valid_genes(expressionfile):
    valid = set()

    with open(expressionfile) as express_in:
        if '.csv' in expressionfile:
            IN = csv.reader(express_in)
        elif '.tsv' in expressionfile:
            IN=csv.reader(express_in, delimiter='\t')
        else:
            print('---ERROR---')
            print('Data filetype not supported')
            exit()

        next(IN)
        for line in IN:
            valid.add(line[0])
    
    return valid