import pandas as pd
import glob
from Bio import SeqIO
from tqdm import tqdm
import os 
import argparse
import time
from functools import partial
import multiprocessing
import ahocorasick

#primers = pd.read_csv("Fw_staggered_primers.tsv", sep='\t', header=None)
#rna = pd.read_csv("gRNAmotifs5.tsv", sep='\t', header=None)
#directory = "/Users/davidchen/Documents/GitHub/Sandor_David/Miniscreen05/Sequence"

def demultiplex(primers, directory, line):
   
    #Ahocorasick automaton
    A = ahocorasick.Automaton()
    for idx, key in enumerate(primers.iloc[:,1]):
        A.add_word(key, (idx, key))
    A.make_automaton()

    primer_dict = dict(zip(primers.iloc[:, 1], primers.iloc[:, 0]))
    
    print(str(line.split("/")[-1].split(".")[0]))
    
    #Iterate through each read in a given fastq file
    with open(line, "r") as handle:
        
        for sequence in tqdm(list(SeqIO.parse(handle, "fastq"))):
        #for sequence in tqdm(SeqIO.parse(handle, "fastq")):
            
            seq_string = str(sequence.seq)
            seq_id = ">" + str(sequence.id)
            
            match_list = []
            for end_index, (insert_order, original_value) in A.iter(seq_string):
                match_list.append(original_value)

            
            if len(match_list) > 0:
                
                with open(directory + "/Split_Fraction/" + line.split("/")[-1].split(".")[0] + "_" + primer_dict[max(match_list , key=len)] + ".fasta", "a+") as file_object:
                    file_object.write(seq_id)
                    file_object.write("\n")
                    file_object.write(seq_string)
                    file_object.write("\n")

            else:
                
                with open(directory + "/Split_Fraction/" + line.split("/")[-1].split(".")[0] + "_unknown.fasta", "a+") as file_object:
                    file_object.write(seq_id)
                    file_object.write("\n")
                    file_object.write(seq_string)
                    file_object.write("\n")
                    
def read_count(primer_sequence_files, rna, directory, fraction):

    read_count_df = pd.DataFrame(data = 0, index =[x.split("/")[-1].split(".")[0]for x in primer_sequence_files], columns = rna.iloc[:, 0])
    rna_dict = dict(zip(rna.iloc[:, 1], rna.iloc[:, 0]))
    
    A = ahocorasick.Automaton()
    for idx, key in enumerate(rna.iloc[:,1]):
        A.add_word(key, (idx, key))
    A.make_automaton()
    
    with open(fraction, "r") as handle:
        
        print(fraction.split("/")[-1].split(".")[0])
        
        #for sequence in tqdm(list(SeqIO.parse(handle, "fasta"))):
        for sequence in tqdm(list(SeqIO.parse(handle, "fasta"))):

            seq_string = str(sequence.seq)
            
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
    parser.add_argument("--primer", "-p", help="Tab separated file of primer sequences to separate based on fraction")
    parser.add_argument("--grna", "-g", help="Tab separated file of gRNA sequences to produce read count file")
    parser.add_argument("--directory", "-d", help="Directory of fastq files for demultiplexing")
    parser.add_argument("--cores", "-c", help="Number of processes (cores) used", default=4)
    # Read arguments from the command line
    args = parser.parse_args()
    
    primers = pd.read_csv(str(args.primer), sep='\t', header=None)
    rna = pd.read_csv(str(args.grna), sep='\t', header=None)
    directory = str(args.directory)
    cores = int(args.cores)
    
    try:
        os.mkdir(directory + "/Split_Fraction")
    except OSError:
        print ("Creation of the directory failed. Check if it already exists!")
    else:
        print ("Successfully created the directory")  
        
    try:
        os.mkdir(directory + "/Split_Fraction_Read_Count")
    except OSError:
        print ("Creation of the directory failed. Check if it already exists!")
    else:
        print ("Successfully created the directory")  
        
    
    sequence_files = sorted([f for f in glob.glob(directory + "/*.fastq")])
    
    demultiplex_iterable = sequence_files
    demultiplex_pool = multiprocessing.Pool(cores)
    demultiplex_func = partial(demultiplex, primers, directory)
    demultiplex_pool.map(demultiplex_func, demultiplex_iterable)
    demultiplex_pool.close()
    demultiplex_pool.join()

    primer_sequence_files = sorted([f for f in glob.glob(directory + "/Split_Fraction/*.fasta")])
    
    read_count_iterable = primer_sequence_files
    read_count_pool = multiprocessing.Pool(cores)
    read_count_func = partial(read_count, primer_sequence_files, rna, directory)
    read_count_pool.map(read_count_func, read_count_iterable)
    read_count_pool.close()
    read_count_pool.join()
    
    read_count_fraction = sorted([f for f in glob.glob(directory + "/Split_Fraction_Read_Count/*.csv")])
    read_count_sum(read_count_fraction, directory)
    
    print('That took {} seconds'.format(time.time() - starttime))


if __name__ ==  '__main__':
    main()       
            
        
        

        
#python read_count.py 
#-p /Users/davidchen/Documents/GitHub/Sandor_David/Miniscreen05/Fw_staggered_primers.tsv 
#-g /Users/davidchen/Documents/GitHub/Sandor_David/Miniscreen05/gRNAmotifs5.tsv 
#-d /Users/davidchen/Documents/GitHub/Sandor_David/Miniscreen05/Sequence  
        
        
#python fast_demultiplex.py -c 8 -p /Users/davidchen/Documents/GitHub/Sandor_David/Miniscreen05/Fw_staggered_primers.tsv -g /Users/davidchen/Documents/GitHub/Sandor_David/Miniscreen05/gRNAmotifs.tsv -d /Users/davidchen/Documents/GitHub/Sandor_David/Miniscreen05/Sequence        
        
