import os, re, time, argparse, sqlite3
import pandas as pd
from multiprocessing import Pool
from Bio.Seq import Seq

def get_cigar_S_len(cigar_string, strand_direction):
    '''
    Calculate soft clipping ('S') lengths from a CIGAR string. In other words,
    number of bases at one end were added sequences not found in the reference.
    
    Input
    -----
    cigar_string: a standard CIGAR (Concise Idiosyncratic Gapped Alignment Report)
    strand_direction: aligned reference strand ('+', '-' or 'NA')
    
    Returns
    -------
    cigar_S_len: integer of soft clipping length.
    
    '''
    
    cigar_S_len = 0 # default
    
    if strand_direction == '+':
        # For '+' strand, find ending soft clipping
        matches = re.findall(r'(\d+)(S{1}$)', cigar_string)
        if matches:
            cigar_S_len = int(matches[0][0])
            
    elif strand_direction == '-':
        # For '-' strand, find beginning soft clipping
        matches = re.findall(r'(^\d+)(S{1})', cigar_string)
        if matches:
            cigar_S_len = int(matches[0][0])
    
    return cigar_S_len

def get_sam_record_tail(sam_record):
    '''
    Extracts SAM info for feature identification.
    
    Input
    -----
    sam_record: a single line from a SAM file
    
    Returns
    -------
    A list of 3 items:
        1. record_id,
        2. soft clippin g length (cigar_S_len),
        3. tail sequence (tail_seq)
            
    '''
    sam_record_tail_params = {}
    sam_record_tail_params['record_id'] = sam_record.split('\t')[0]
    cigar_string = sam_record.split('\t')[5]
    sam_seq = sam_record.split('\t')[9]
    
    if sam_record.split('\t')[1] == '0':
        sam_record_tail_params['strand'] = '+'
        cigar_S_len = get_cigar_S_len(cigar_string, '+')
        sam_record_tail_params['cigar_S_len'] = cigar_S_len
        tail_seq = sam_seq[-cigar_S_len:]
        sam_record_tail_params['tail_seq'] = tail_seq
        
    elif sam_record.split('\t')[1] == '16':
        sam_record_tail_params['strand'] = '-'
        cigar_S_len = get_cigar_S_len(cigar_string, '-')
        sam_record_tail_params['cigar_S_len'] = cigar_S_len
        tail_seq = sam_seq[:cigar_S_len]
#         sam_record_tail_params['tail_seq'] = tail_seq
        sam_record_tail_params['tail_seq'] = str(Seq(tail_seq).reverse_complement())
        
    else:
        sam_record_tail_params['strand'] = 'NA'
        sam_record_tail_params['cigar_S_len'] = 0
        sam_record_tail_params['tail_seq'] = ''
    
    result = [sam_record_tail_params['record_id'],
#               sam_record_tail_params['strand'],
              sam_record_tail_params['cigar_S_len'],
              sam_record_tail_params['tail_seq']
             ]
    
    return result

def get_sam_tails(SAMfile, n_processes=4, output_folder=None):
    '''
    Identify the tail (unmapped part) of each mapped sequence in the SAMfile.
    
    Input
    -----
    SAMfile:
        path to the .sam alignment file
    n_processes:
        number of processes for parallel processing; default 4
    output_folder:
        folder to store output files; default is the parental folder of SAMfile
    
    Output
    ------
    A csv file with matching file name of the SAMfile in the output folder.
    
    The csv file has 3 columns:
        record_id, cigar_S_len and tail_seq.
        
    '''
    # Record starting time for benchmarking time cost
    start = time.time()
    
    # Generate a list of paired parameters to parallel process
    lines = []
    with open(SAMfile) as f:
        for line in f:
            if not line.startswith('@'):
                lines.append(line)

    # Use multpile worker processes to process all entries
    with Pool(processes=n_processes) as pool:
        output = pool.map(get_sam_record_tail, lines)

    # Generate the output csv file name to match the SAM file name
    if output_folder is None:
        output_folder = os.path.dirname(SAMfile)
    SAM_base_name = os.path.basename(SAMfile)[:-4]
    output_csv_file = os.path.join(output_folder, SAM_base_name+'_tails.csv')
    
    # Store the output
    df_output = pd.DataFrame(output, columns=['record_id', 'cigar_S_len', 'tail_seq'])
    df_output.to_csv(output_csv_file, index=False)

    end = time.time()
    print('Time cost processing', len(lines), 'SAM records is:', end - start, 'seconds.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("SAMfile",
                        help="path to the .sam alignment file")
    parser.add_argument("n_processes", nargs='?',
                        help="optional; the number of processes to use",
                        default=4)
    
    args = parser.parse_args()
    get_sam_tails(args.SAMfile, args.n_processes)