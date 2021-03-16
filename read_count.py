import pandas as pd
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm
import os 
import argparse
import time
from functools import partial
import multiprocessing
import ahocorasick

#rna = pd.read_csv("gRNAlist.tsv", sep='\t', header=None)
#directory = "/Users/davidchen/Documents/GitHub/Sandor_David/ARenh"


def read_count(primer_sequence_files, rna, directory, trim_start, trim_end, fraction):

    indices = sorted([f for f in glob.glob(directory + "/*.fastq")])
    
    read_count_df = pd.DataFrame(data = 0, index =[x.split("/")[-1].split(".")[0]for x in indices], columns = rna.iloc[:, 0])
    rna_dict = dict(zip(rna.iloc[:, 1], rna.iloc[:, 0]))
    
    A = ahocorasick.Automaton()
    for idx, key in enumerate(rna.iloc[:,1]):
        A.add_word(key, (idx, key))
    A.make_automaton()
    
    with open(fraction, "r") as handle:
        
        print(fraction.split("/")[-1].split(".")[0])
        
        #for sequence in tqdm(list(SeqIO.parse(handle, "fasta"))):
        for sequence in tqdm(list(SeqIO.parse(handle, "fastq"))):

            seq_string = str(sequence.seq)[trim_start:trim_end]
            
            for end_index, (insert_order, original_value) in A.iter(seq_string):
                
                #assertion
                start_index = end_index - len(original_value) + 1
                assert seq_string[start_index:start_index + len(original_value)] == original_value
                
                read_count_df.loc[fraction.split("/")[-1].split(".")[0], rna_dict[original_value]] += 1
            
    read_count_df.to_csv(directory + "/Split_Fraction_Read_Count/" + fraction.split("/")[-1].split(".")[0] + "_read_counts.csv", index=True)
    

def read_count_sum(read_count_fraction, directory):
    
    df = pd.read_csv(read_count_fraction[0], index_col=0)
    
    for fraction_reads in range(1, len(read_count_fraction)):
        df += pd.read_csv(read_count_fraction[fraction_reads], index_col=0)
        
    df.to_csv(directory + "/read_counts.csv", index=True)
    
def main(): 
    
    starttime= time.time()
    
    # Initiate the parser
    parser = argparse.ArgumentParser(description="Demultiplex Parameters")
    
    # Add long and short argument
    parser.add_argument("--grna", "-g", help="Tab separated file of gRNA sequences to produce read count file")
    parser.add_argument("--directory", "-d", help="Directory of fastq files for demultiplexing")
    parser.add_argument("--cores", "-c", help="Number of processes (cores) used", default=4)
    parser.add_argument("--paired", "-p", help="True = Paired; False = Single", default=False)
    parser.add_argument("--forward_id", "-a", help="ID string in file name of forward reads", default=None)
    parser.add_argument("--reverse_id", "-b", help="ID string in file name of reverse reads", default=None)
    parser.add_argument("--trim_start", "-s", help="Trim read to start at this index (first index = 0)", default=None)
    parser.add_argument("--trim_end", "-e", help="Trim read to end at this index", default=None)
    # Read arguments from the command line
    args = parser.parse_args()
    
    forward_rna = pd.read_csv(str(args.grna), sep='\t', header=None)
    reverse_rna = forward_rna.copy()
    reverse_rna.iloc[:,1] = [str(Seq(x).reverse_complement()) for x in reverse_rna.iloc[:,1]]
    
    directory = str(args.directory)
    cores = int(args.cores)
    paired = str(args.paired)
    
    forward_id = str(args.forward_id)
    reverse_id = str(args.reverse_id)
    trim_start = int(args.trim_start)
    trim_end = int(args.trim_end) + 1

    try:
        os.mkdir(directory + "/Split_Fraction_Read_Count")
    except OSError:
        print ("Creation of the directory failed. Check if it already exists!")
    else:
        print ("Successfully created the directory")  
        

    primer_sequence_files = sorted([f for f in glob.glob(directory + "/*.fastq")])
    
    forward_read_files = [x for x in primer_sequence_files if forward_id in x]
    reverse_read_files = [x for x in primer_sequence_files if reverse_id in x]
    
    if paired == "True":
        
        print("Paired End")
        read_count_iterable = forward_read_files
        read_count_pool = multiprocessing.Pool(cores)
        read_count_func = partial(read_count, forward_read_files, forward_rna, directory, trim_start, trim_end)
        read_count_pool.map(read_count_func, read_count_iterable)
        read_count_pool.close()
        read_count_pool.join()
        
        read_count_iterable = reverse_read_files
        read_count_pool = multiprocessing.Pool(cores)
        read_count_func = partial(read_count, reverse_read_files, reverse_rna, directory, trim_start, trim_end)
        read_count_pool.map(read_count_func, read_count_iterable)
        read_count_pool.close()
        read_count_pool.join()
        
    if paired == "False":
        print("Single End")
        read_count_iterable = primer_sequence_files
        read_count_pool = multiprocessing.Pool(cores)
        read_count_func = partial(read_count, primer_sequence_files, forward_rna, directory, trim_start, trim_end)
        read_count_pool.map(read_count_func, read_count_iterable)
        read_count_pool.close()
        read_count_pool.join()
        
        
    read_count_fraction = sorted([f for f in glob.glob(directory + "/Split_Fraction_Read_Count/*.csv")])
    read_count_sum(read_count_fraction, directory)
    
    print('That took {} seconds'.format(time.time() - starttime))


if __name__ ==  '__main__':
    main()       
    
        
        
#python read_count.py -c 4 -g /Users/davidchen/Documents/GitHub/Sandor_David/ARenh/gRNAlist.tsv -d /Users/davidchen/Documents/GitHub/Sandor_David/ARenh -p True -a R1 -b R2      
        
