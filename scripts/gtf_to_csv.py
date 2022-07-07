import numpy as np
import pandas as pd
import os
import argparse

def get_gene_bed(GTFfile, gene_bed_file=''):
    '''
    Inherited from gtftools.py. Removed chromosome filtering.
    
    Input
    -----
    GTFfile: path to the gtf file containing all gene annotations
    gene_bed_file: file path to store the output bed file
    
    Output
    ------
    If gene_bed_file is provided, writes the gene annotation in bed format
     
    Returns
    -------
    genebed: a dictionary of the gene_bed info
    
    '''

    genebed={}
    f = open(GTFfile)
    for line in f:
        table = line.split('\t')
        if line[0] != '#':  # skip comment line
            if table[2] == "gene":
                ensid  = line.split('gene_id')[1].split('"')[1]
                symbol = line.split('gene_name')[1].split('"')[1]
                record = (table[0],int(table[3])-1,int(table[4]),table[6],ensid,symbol)
                if table[0] in genebed:
                    genebed[table[0]].append(record)
                else:
                    genebed[table[0]] = [record]
    f.close()

    # print to file if required
    if len(gene_bed_file) > 1:
        f=open(gene_bed_file,'w')
        for ichrom in genebed.keys():
            ichrom_gene = genebed[ichrom]
            for igene in ichrom_gene:
                f.write('\t'.join([igene[0],str(igene[1]),str(igene[2]),igene[3],igene[4],igene[5],])+'\n')
        f.close()

    return genebed

def bedmerge(featureRange):
    # featureRange: a list of ranges in bed format, such as [(1,1000,2000,+),(1,2200,3000,-)]
    featureRange.sort(key = lambda x: (x[0],x[1]))
    merged=[]
    nRange=len(featureRange)
    if nRange == 1:
        merged=featureRange
    elif nRange == 2:
        imerge=neighbor_merge(featureRange[0],featureRange[1])
        for each in imerge:
            merged.append(each)
    else:
        i = 2
        imerge=neighbor_merge(featureRange[0],featureRange[1])
        n_imerge=len(imerge)
        while n_imerge > 0:
            if n_imerge == 2:
                merged.append(imerge[0])

            imerge=neighbor_merge(imerge[n_imerge-1],featureRange[i])
            n_imerge=len(imerge)
            if i == nRange-1:
                for each in imerge:
                    merged.append(each)
                n_imerge = -1
            i+=1

    return merged

def neighbor_merge(range1, range2):
    if range2[1]<=range1[2]:
        merged =[(range1[0],range1[1],max(range1[2],range2[2]),range1[3])]
    else:
        merged=[range1,range2]
    return merged

def merge_exon(GTFfile, merged_exon_file=''):
    # record exon coordination
    
    exon={}
    f=open(GTFfile)
    for line in f:
        table=line.split('\t')
        if line[0] != '#':  # skip comment line
            if table[2] == 'exon':
                gene=line.split('gene_id')[1].split('"')[1]
                iexon=(table[0],int(table[3])-1,int(table[4]),table[6])  # gtf to bed coordination
                if gene in exon:
                    exon[gene].append(iexon)
                else:
                    exon[gene]=[iexon]
    f.close()
    # merge all exons of each gene
    merged_exon={}
    gene_length={}
    for gene in exon:
        merged=bedmerge(exon[gene])
        merged_exon[gene]=merged

        # calculate merged gene length(sum of non-overlapping exons)
        length=0
        for each in merged:
            length+=each[2]-each[1]
        gene_length[gene]=length

    # assign merged exons and gene lengths
    merged = {}
    merged['merged_exon']=merged_exon
    merged['merged_gene_length']=gene_length


    # print merged exons if required
    if len(merged_exon_file) > 1:
        m=open(merged_exon_file,'w')
        for gene in merged_exon:
            imerged=merged_exon[gene]
            for each in imerged:
                joined="\t".join([each[0],str(each[1]),str(each[2]),gene,'0',each[3]])+"\n"
                m.write(joined)
        m.close()

    return merged

def get_isoform_length(GTFfile, isoformlength_file=''):
    # calculate length of isoforms
        
    f=open(GTFfile)
    isolength = {}
    gene2iso= {}
    iso2gene= {}
    for line in f:
        table=line.split('\t')
        if line[0] != '#':  # skip comment line
            if table[2] == 'exon':
                gene= line.split('gene_id')[1].split('"')[1]  # gene ID
                tcx = line.split('transcript_id')[1].split('"')[1]  # tcx ID
                exon_length = int(table[4])-int(table[3])+1

                # isoform length calculate
                if tcx in isolength:
                    isolength[tcx] = isolength[tcx] + exon_length
                else:
                    isolength[tcx] =  exon_length

                # record isoforms of genes
                if gene in gene2iso:
                    gene2iso[gene].append(tcx)
                else:
                    gene2iso[gene] = [tcx]
                iso2gene[tcx] = gene

    if len(isoformlength_file) > 1:
        f=open(isoformlength_file,'w')
        f.write('isoform\tgene\tlength\n')
        for thisiso in iso2gene:
            out=thisiso+'\t'+iso2gene[thisiso]+'\t'+str(isolength[thisiso])+'\n'
            f.write(out)
        f.close()

    # return
    ret={}
    ret['gene2iso']  = gene2iso
    ret['isoform_length'] = isolength
    return(ret)

def get_gene_length(GTFfile, genelength_file=''):
    # record exon coordination

    # calculate length of merged exons
    merged_gene_length = merge_exon(GTFfile)['merged_gene_length']

    # calculate length of each transcrpt isoform
    ret_isoform_length = get_isoform_length(GTFfile)
    gene2iso =ret_isoform_length['gene2iso']
    isolength=ret_isoform_length['isoform_length']

    # calculate gene length as mean, median or max of isoforms
    gene_length = {}
    for igene in gene2iso:
        isoforms = gene2iso[igene]
        length_set = []
        for iisoform in isoforms:
            length_set.append(isolength[iisoform])
        gene_length[igene] = [int(np.mean(length_set)),int(np.median(length_set)),max(length_set),merged_gene_length[igene]]

    # print
    if len(genelength_file) > 1:
        f = open(genelength_file,'w')
        f.write('gene_id\tmean\tmedian\tlongest_isoform\tmerged\n')
        for igene in gene_length:
            tmp =  gene_length[igene]
            tmplst = [igene,str(tmp[0]),str(tmp[1]),str(tmp[2]),str(tmp[3])]
            f.write('\t'.join(tmplst)+'\n')
        f.close()
        
def gtf_to_csv(GTFfile):
    '''
    Input
    -----
    GTFfile: path to a .gtf gene annotation file
    
    Output
    ------
    Generate 3 files in the same folder with matching base filenames
        1. a *_genes.bed file with gene position info
        2. a *_gene_length.csv gene length file (4 variants of gene lengths provided)
        3. a *_genes.csv file combining gene position and length info
    '''
    parent_folder = os.path.dirname(GTFfile)
    base_name = os.path.basename(GTFfile)[:-4]
    gene_bed_file = os.path.join(parent_folder, base_name+'_genes.bed')
    gene_length_file = os.path.join(parent_folder, base_name+'_gene_length.csv')
    genes_csv_file = os.path.join(parent_folder, base_name+'_genes.csv')
    
    if os.path.isfile(genes_csv_file):
        print('A matching genes.csv file is already present.\nRemove it if you want to regenerate from GTF.')
        return
    
    if not os.path.isfile(gene_bed_file):
        print('No matching gene bed file found.\nGenerating gene bed file...')
        get_gene_bed(GTFfile, gene_bed_file)
    else:
        print('Found matching gene bed file.\nWill use it for gene bed info.')

    if not os.path.isfile(gene_length_file):
        print('No matching gene length file found.\nGenerating gene length file...')
        get_gene_length(GTFfile, gene_length_file)
    else:
        print('Found matching gene length file.\nWill use it for gene length info.')

    df_gene_bed = pd.read_csv(gene_bed_file, sep='\t', header=None)
    df_gene_bed.columns = ['chr', 'start', 'end', 'strand', 'gene_id', 'gene_name']

    df_gene_length = pd.read_csv(gene_length_file, sep='\t')
    assert 'gene_id' in df_gene_length.columns

    df_genes = pd.merge(df_gene_bed, df_gene_length, on='gene_id')
    # Check expected columns are present in the data frame
    for i in ['chr', 'start', 'end', 'strand', 'gene_id', 'gene_name',
              'mean', 'median', 'longest_isoform', 'merged']:
        assert i in df_genes.columns
    df_genes.to_csv(genes_csv_file, index=False)
    
    print('Successfully generated genes.csv file from GTF.')
    
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("GTFfile",
                        help="path to the .gtf gene annotation file")
    args = parser.parse_args()
    gtf_to_csv(args.GTFfile)