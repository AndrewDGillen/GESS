import subprocess as sub
import os
from collections import defaultdict
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy

def process_genelist(genelist, seqfile):
    allgenes = []

    with open(genelist) as genesIN:
        for line in genesIN:
            allgenes.append(line.strip())

    geneseq = ''
    genename = None

    seq_by_gene = {}
    maxlen_by_gene = defaultdict(lambda:0)

    with open(seqfile) as seqsIN:
        for line in seqsIN:
            if line[0] == '>':
                
                if genename in allgenes:
                    if len(geneseq) > maxlen_by_gene[genename]:
                        maxlen_by_gene[genename] = len(geneseq)
                        seq_by_gene[genename] = geneseq

                geneseq = ''

                seqname = line.split('name=')[1].split(';')[0]
                genename = seqname[:-3]

            else:
                geneseq += line.strip()

        if genename in allgenes:
            if len(geneseq) > maxlen_by_gene[genename]:
                maxlen_by_gene[genename] = len(geneseq)
                seq_by_gene[genename] = geneseq

    try:
        assert len(allgenes) == len(seq_by_gene)
    except:
        print("ERROR: Some genes of interest don't have associated sequences in the provided file.")
        print(len(allgenes))
        print(len(seq_by_gene))#
        exit()

    return allgenes, seq_by_gene

def pairwise_align_seqs(seq1, seq2, comparegene):

    with open('seq1.fa', 'w+') as seq1file:
        seq1file.write(seq1)
    with open('seq2.fa', 'w+') as seq2file:
        seq2file.write(seq2)

    water_output = sub.run(('needle', '-asequence', 'seq1.fa', '-bsequence','seq2.fa', '-sprotein1', '-sprotein2', '-gapopen', '10.0', '-gapextend', '0.5', '-outfile', 'stdout'), capture_output=True)
    os.remove('seq1.fa')
    os.remove('seq2.fa')

    water_output = water_output.stdout.decode("utf-8")

    identity = water_output.split('Identity:')[1].split('(')[1].split('%)\n')[0]
    similarity = water_output.split('Similarity:')[1].split('(')[1].split('%)\n')[0]

    return float(identity), float(similarity)

def handle_gene_matches(i_seq, genesequences, all_genes):

    identity_vector = []
    similarity_vector = []

    for comparegene in all_genes:
        compareseq = genesequences[comparegene]
        identity, similarity = pairwise_align_seqs(i_seq, compareseq, comparegene)

        identity_vector.append(identity)
        similarity_vector.append(similarity)
    
    return identity_vector,similarity_vector

def main(genelist, seqfile):
    all_genes, genesequences = process_genelist(genelist, seqfile)

    all_intergene_ids = []
    all_intergene_sims = []

    for genename in all_genes:
        i_seq = genesequences[genename]

        ids, sims = handle_gene_matches(i_seq, genesequences, all_genes)

        all_intergene_ids.append(ids)
        all_intergene_sims.append(sims)

    id_array = np.array(all_intergene_ids)
    sim_array = np.array(all_intergene_sims)

    for harvest, filename in [(id_array, 'IDMatrix.svg'), (sim_array, 'SimMatrix.svg')]:


        ordered_columns = scipy.cluster.hierarchy.linkage(np.transpose(harvest), method='weighted', metric='euclidean', optimal_ordering=True)
        ordered_rows = scipy.cluster.hierarchy.linkage(harvest,method='weighted', metric='euclidean', optimal_ordering=True)
    
        #Generates a pandas datafame containing the gess scores, labelled columns and rows as appropriate
        df = pd.DataFrame(harvest, columns = all_genes)
        df.index= all_genes

        print(df)
        #Generates a Seaborn Clustermap from the supplied data. Parameters:
            #Colour range 0-100 (as GESS are percentages)
            #Labels on each row and column
            #Clustering for both rows and columns based on the user-selected clustering algorith, 
            #   by default we use the Weighted clustering algorithm/ WPGMA (Weighted Pair Group Method with Arithmetic Mean) 
            #Uses pre-generated linkage matrices to preserve optimal (ie minimal) euclidean distance between neighbors
        
        #Note that this clustering CAN use fastcluster, but doesn't NEED to
        sns.clustermap(df, vmin=0, vmax=100, annot=False, method='weighted', xticklabels = 1, yticklabels=1, row_cluster=True, col_cluster=True, row_linkage = ordered_rows, col_linkage=ordered_columns)
        plt.savefig(filename)


if __name__ == '__main__':
    main('/home/ag_jdowlab/Desktop/FlyGene/ALL_vha_genes.txt', '/home/ag_jdowlab/Desktop/FlyGene/VHA_seqs_proteins.fasta')
