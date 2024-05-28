import argparse

def get_gess_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-qg', '--querygene', type = str, default='', help = 'Query Gene')
    parser.add_argument('-tg', '--targetgene', type = str, default='', help = 'Target Gene')
    parser.add_argument('-qd', '--querydata', type=str, default='', help ='Enrichment data corresponding to the Query genes')
    parser.add_argument('-td', '--targetdata', type=str, default='', help ='Enrichment data corresponding to the Target genes')
    parser.add_argument('-qs', '--queryspecies', type=str, default = '', help ='The species from which query genes are derived. This only needs to be defined if not a default FlyGene species')
    parser.add_argument('-ts', '--targetspecies', type=str, default = '', help ='The species from which target genes are derived. This only needs to be defined if not a default FlyGene species')
    parser.add_argument('-a', '--annotations', type=str, default='', help = 'Tissue annotation file. Required if using a flat expression matrix')
    parser.add_argument('-test', '--testing', type=bool, default=False, help = 'Boolean switch to activate GESS testing suite.')
    parser.add_argument('-v', '--verbose', type = bool, default = False, help = 'Boolean switch activating VERBOSE mode')
    parser.add_argument('-u', '--unacceptable', type=str, default='', help="file containing uninterpretable annotations. Line seperated text file (see example)")
    parser.add_argument('-anno', '--targetannos', nargs='*', default = [], help = "Selected annotation levels to use in the GESS calculation")
    parser.add_argument('-mode', '--analysis_mode', type=str, choices = [None, 'bulk', 'prevalence','expression'], default = None, help = 'Mode of h5ad analysis - either "prevalence" (the % of reads from a context with detected expression UMI>1) or "expression" (the average expression within a given context)')
    parser.add_argument('-umi', '--umithreshold', type = int, default = 1, help = 'Sensitivity threshold for H5 parsing in expression mode. Any cell expressing a gene >= this threshold will count as an "Expressing" cell')
    parser.add_argument('-s', '--save', type = str, default='', help = 'Filename under which to save the produced matrix. Please only use image file formats. Leave blank if you do not wish to save')
    parser.add_argument('-c', '--cluster_method', type = str, choices = ['single', 'complete', 'median', 'average', 'weighted', 'centroid', 'ward'], default = 'weighted', help = 'Clustering method for heatmap/dendogram generation')
    parser.add_argument('-l', '--label', type=bool, default=False, help = 'Boolean switch for activating labelling of the heatmap. This gets very busy for big heatmaps!')
    parser.add_argument('-iq', '--input_query', type=str, default = '', help = 'Chosen genes from the query species you wish to include within the matrix')
    parser.add_argument('-it', '--input_target', type=str, default='', help = 'Chosen genes from the target species you wish to include within the matrix')
    args = parser.parse_args()

    return args