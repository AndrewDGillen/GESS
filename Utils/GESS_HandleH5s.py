import pandas as pd
import scanpy as sp
import gc
from collections import Counter

#Returns all groupable annotations from the H5 file - these are the possible annotation levels for GESS
#Note that annotation levels with only one category cannot be used for GESS.
#Levels with only one category are returned in "Other_annos" to inform error printing
def get_annotation_levels(anndata):

    annotations = list(anndata.obs.columns)

    usable_annos = []
    other_annos = {}

    for j in annotations:

        if len(anndata.obs[j].unique()) > 1:
            usable_annos.append(j)
        else:
            other_annos[j] =  list(anndata.obs[j].unique())[0]

    return usable_annos, other_annos

#finds "Expressing Cells" from a single gene's expression data in a whole dataset.
#These are summed to provide a single count.
def get_expressing_cells(subarray, threshold):
    
    return sum([1 for cell in subarray if cell >= threshold ])

#The major function for handling H5s.
#Given a h5 file, a set of annotations, a gene of interest, an analysis mode (expression or prevalence, see below), and optionally a UMI threshold, returns a sorted data vector
#This data vector is a dictionary, containing a sorted list of annotations at each requested annotation level (higher/more expression to lower/least)

#ANALYSIS MODES:
    #Expression - Ranks annotations based on average expression of the gene of interest
    #Prevalence - Ranks annotations based on the proportion of cells which express a gene of interest >= the UMI threshold

#The UMI threshold (default 1) is the level of expression required to define a cell as expressing a gene
def H5_GESS(desired_annotations, designated_file, gene_of_interest, analysis_mode, umithresh = 1):

    #Uses ScanPy to load a H5 file
    adata, frame_settings = prep_h5(designated_file)

    #finds usable annotations from the supplied h5 file
    usable_annotations, other_annos = get_annotation_levels(adata) 
    
    #Gets a list of all genes within the H5 files
    all_genes = frame_settings.var_names.to_list()

    #Error handling for cases where the gene of interest is not found in the supplied H5 file
    if gene_of_interest not in all_genes:
        print(f'--GESS ERROR--\nGene of Interest {gene_of_interest} not in the supplied h5 file. Please check naming conventions are consistent!')
        exit()
    
    #Deleted to optimise memory usage
    del all_genes 

    #Prepares output results
    data_by_annotation_level = {}

    #This is the most efficient (time AND memory) way of filtering large adata objects to retrieve only the column associated with the gene of interest
    adata_formatted = frame_settings[:, gene_of_interest].X.toarray()

    #Each specified annotation level is handled seperately, providing the multi-layered data required for GESS
    for annotation in desired_annotations:

        #Error handling for cases where the user wants to use an annotation which doesn't exist
        if annotation not in usable_annotations:

            if annotation in other_annos:
                print(f'--GESS ERROR--\nDesired annotation "{annotation}" only has one category ({other_annos[annotation]}) and so cannot be used for GESS. Please choose other annotations.')

            else:
                print(f'--GESS ERROR--\nDesired annotation {annotation} not in the grouping attributes of the supplied h5 file')

            print(f'Appropriate attributes:')

            for j in usable_annotations:
                print(f'\t{j}')

            exit()

        #Gets all possible contexts associated with the desired annotation
        cell_annotations = adata.obs[annotation].to_list()
        unique_annos = set(cell_annotations)

        #Moves the adata object to a pandas dataframe for ease of processing
        formatted = pd.DataFrame(adata_formatted, index=cell_annotations, columns = ['GOI'])
        formatted.index.name = 'Annotation'

        #If the analysis mode is expression, we find the Mean expression for each annotation type
        if analysis_mode == 'expression':
            annotation_means = formatted.groupby('Annotation').mean().to_dict()['GOI']
            gene_data_vector = [(x, annotation_means[x]) for x in annotation_means]

            for i_annotation in unique_annos:
                if i_annotation not in annotation_means:
                    gene_data_vector.append((i_annotation, 0))

        #If the analysis mode is prevalence, we find the proportion of cells of each annotation type which express a gene over the UMI threshold
        elif analysis_mode == "prevalence":

            #Gathers the total number of cells belonging to each annotation type
            totals = Counter(cell_annotations)
            
            #Filters the DF to keep only rows in which the gene is expresed > UMI threshold
            formatted = formatted[formatted.GOI >= umithresh]

            #Counts the number of cells remaining in each annotation type
            annotation_counts = formatted.groupby('Annotation').count().to_dict()['GOI']

            #Returns the proportion of expressing cells to the data vector
            gene_data_vector = [(x, (annotation_counts[x]/totals[x])) for x in annotation_counts]
            
            # print(gene_data_vector)

            for i_annotation in unique_annos:
                if i_annotation not in annotation_counts:
                    gene_data_vector.append((i_annotation, 0))

        # print(gene_data_vector)
        data_by_annotation_level[annotation]=sorted(gene_data_vector, key= lambda x:x[1], reverse=True)

    del adata
    gc.collect()
    
    return data_by_annotation_level

#A utility feature allowing users to interrogate H5AD files
def describe_file(designated_file):
    adata, frame_settings = prep_h5(designated_file)

    usable_annotations, other_annotations = get_annotation_levels(adata) 

    print("Number of cells:", adata.n_obs)
    print("Number of genes:", len(frame_settings.var_names))

    #User-friendly printing of AnnData attributes to help in selection
    print('Appropriate attributes:')
    for j in usable_annotations:

        unique_annos = [str(i) for i in adata.obs[j].unique()] 
        short_annos = unique_annos[:min(len(unique_annos), 5)]

        if len(unique_annos) > 1:
            print(f'--{j}')
            print(f'For example: {"; ".join(short_annos)}')
            print()

    #We also let the users know which Annotations are inappropriate for GESS and why.
    if len(other_annotations) > 0:

        print('Some annotation levels have only one category - GESS cannot utilise these:')
        for level in other_annotations:
            print(f'{level}\t Only Category:   {other_annotations[level]}')

def test():
    #testfile = '/home/ag_jdowlab/r_fca_biohub_malpighian_tubule_10x.h5ad'
    #testfile = '/home/ag_jdowlab/r_fca_biohub_all_wo_blood_10x.h5ad'
    testfile = '/media/ag_jdowlab/Expansion/HandlingLFC/UpdatedH5s/MsIT_Updated.h5ad'
    testanno = ["annotation", "annotation_broad"]
    testgenes = 'LOC115452582'
    analysis_mode = 'prevalence'

    describe_file(testfile)
    #H5_GESS(testanno, testfile, testgenes, analysis_mode)

#This function is designed to try and handle a more broad range of H5 files
#We aim to find where expression data and gene names are stored.
def prep_h5(infile):

    adata = sp.read_h5ad(infile)

    frame_setting = adata

    if adata.raw != None:
        frame_setting = adata.raw

    return adata, frame_setting

def get_h5_genes(h5file):

    adata, frame_setting = prep_h5(h5file)

    #Gets a list of all genes within the H5 files
    all_genes = frame_setting.var_names.to_list()

    return all_genes

if __name__ == '__main__':
    test()