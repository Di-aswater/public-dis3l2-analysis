import os, re, time, argparse
import pandas as pd
from gtf_to_csv import gtf_to_csv
from multiprocessing import Pool

def get_cigar_len(cigar_string):
    '''
    Calculate seq lengths from a CIGAR string
    
    Input
    -----
    cigar_string: a standard CIGAR (Concise Idiosyncratic Gapped Alignment Report)
    
    Returns
    -------
    cigar_len: a dictionary of 2 items:
               temp_len, the consumed length on the alignment template.
               aln_seq_len, the length of sequence aligned to the template.
               Note: not counting insertions and softclipping.
    '''
    # use regular expression to parse a given cigar string
    matches = re.findall(r'(\d+)([A-Z]{1})', cigar_string)
    
    template_len = 0
    aln_seq_len = 0
    for m in matches:
        # Validate CIGAR string
        assert m[1] in ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']
        
        # Add up consumed template length
        if m[1] in ['M', 'D', 'N', '=', 'X']:
            template_len += int(m[0])
        
        # Add up aligned seq length (not counting soft clipping)
        if m[1] in ['M', '=', 'X']:
            aln_seq_len += int(m[0])
    
    # Compose a dictionary to return
    cigar_len = {'template_len': template_len, 'aln_seq_len': aln_seq_len}
    
    return cigar_len

def get_sam_record_params(sam_record):
    '''
    Extracts SAM info for feature identification.
    
    Input
    -----
    sam_record: a single line from a SAM file
    
    Returns
    -------
    sam_record_params: a dictionary of SAM info for feature identification, including
                        record_id,
                        alignment flag (strand),
                        chromosome (aln_chr),
                        starting chr position (aln_pos; included),
                        ending chr position (aln_pos_end; included),
                        and aligned length (aln_seq_len).
    '''
    sam_record_params = {}
    sam_record_params['record_id'] = sam_record.split('\t')[0]
    if sam_record.split('\t')[1] == '0':
        sam_record_params['strand'] = '+'
    elif sam_record.split('\t')[1] == '16':
        sam_record_params['strand'] = '-'
    else:
        sam_record_params['strand'] = 'NA'
    sam_record_params['aln_chr'] = sam_record.split('\t')[2]
    sam_record_params['aln_pos'] = int(sam_record.split('\t')[3]) # 1-based, included
    record_cigar = sam_record.split('\t')[5]
    cigar_len = get_cigar_len(record_cigar)
    sam_record_params['aln_seq_len'] = cigar_len['aln_seq_len']
    sam_record_params['aln_pos_end'] = sam_record_params['aln_pos'] + cigar_len['template_len'] - 1

    return sam_record_params

def id_gene(df_genes, sam_record):
    '''
    Identify the gene  SAM info for feature identification.
    
    Input
    -----
    df_genes: a pandas data frame containing gene position and length info
    sam_record: a single line from a SAM file
    
    Returns
    -------
    a list of 3 items:
        1. record_id from the sam_record
        2. gene_id that was uniquely matched; 'NA' if no match or multiple matches
        3. gene_name that was uniquely matched; 'NA' if no match or multiple matches
    '''

    # Get parameters required for matching gene id from the sam record
    params = get_sam_record_params(sam_record)
    
    # Narrow down the data frame for possible matches
    df_search = df_genes[(df_genes.chr == params['aln_chr']) &
                         (df_genes.strand == params['strand']) &
                         (df_genes.start <= params['aln_pos']) &
                         (df_genes.end >= params['aln_pos_end']) &
                         (df_genes.longest_isoform >= params['aln_seq_len'])]
    
    if len(df_search) == 1:
        return [params['record_id'], df_search.gene_id.tolist()[0], df_search.gene_name.tolist()[0]]
    else:
        return [params['record_id'], 'NA', 'NA']

def f_id_gene(x):
    '''
    A wrapper function to accept a list of parameters for parallel processing
    
    Input
    -----
    x[0]: a pandas data frame containing gene position and length info
    x[1]: a single line of SAM record
    
    Returns
    -------
    a list of 3 items:
        1. record_id from the sam_record
        2. gene_id that was uniquely matched; 'NA' if no match or multiple matches
        3. gene_name that was uniquely matched; 'NA' if no match or multiple matches
    '''
    return id_gene(x[0], x[1])

def id_genes_sam(SAMfile, GTFfile, n_processes=4):
    '''
    Identify genes of each mapped sequence in the SAMfile using
    GTFfile as the reference annotation.
    
    Input
    -----
    SAMfile: path to the .sam alignment file
    GTFfile: path to the .gtf gene annotation file
    n_processes: number of processes for parallel processing
    
    Output
    ------
    A csv file with matching file name in the same folder of the SAMfile
    with 3 columns: record_id, gene_id and gene_name
    '''
    
    # Record starting time for benchmarking time cost
    start = time.time()
    
    # Generate genes.csv if not already exists
    GTF_parent_folder = os.path.dirname(GTFfile)
    GTF_base_name = os.path.basename(GTFfile)[:-4]
    genes_csv_file = os.path.join(GTF_parent_folder, GTF_base_name+'_genes.csv')
    if not os.path.isfile(genes_csv_file):
        gtf_to_csv(GTFfile)

    # Read in the data frame from genes.csv
    df_genes = pd.read_csv(genes_csv_file)

    # Generate a list of paired parameters to parallel process
    lines = []
    with open(SAMfile) as f:
        for line in f:
            if not line.startswith('@'):
                lines.append(line)
    param_list = [[df_genes, line_] for line_ in lines]
    
    # Use multpile worker processes to process all entries
    with Pool(processes=n_processes) as pool:
        output = pool.map(f_id_gene, param_list)

    # Generate the output csv file name to match the SAM file name
    SAM_parent_folder = os.path.dirname(SAMfile)
    SAM_base_name = os.path.basename(SAMfile)[:-4]
    output_csv_file = os.path.join(SAM_parent_folder, SAM_base_name+'_out.csv')
    
    # Store the output
    df_output = pd.DataFrame(output, columns=['record_id', 'gene_id', 'gene_name'])
    df_output.to_csv(output_csv_file)

    end = time.time()
    print('Time cost processing', len(lines), 'SAM records is:', end - start, 'seconds.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("SAMfile",
                        help="path to the .sam alignment file")
    parser.add_argument("GTFfile",
                        help="path to the .gtf gene annotation file")
    parser.add_argument("n_processes", nargs='?',
                        help="optional; the number of processes to use",
                        default=4)
    args = parser.parse_args()
    id_genes_sam(args.SAMfile, args.GTFfile, args.n_processes)