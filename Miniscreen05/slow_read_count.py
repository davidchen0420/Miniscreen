import pandas as pd
import glob
from Bio import SeqIO
import regex
from tqdm import tqdm
import os 
import argparse
import csv


primers = pd.read_csv("Fw_staggered_primers.tsv", sep='\t', header=None)
rna = pd.read_csv("gRNAmotifs5.tsv", sep='\t', header=None)
directory = "/Users/davidchen/Documents/GitHub/Sandor_David/Miniscreen05/Sequence"
n_mismatch = 1

def demultiplex(primers, directory, m_mismatch = 0):

    #Iterate through all fastq read files in directory
    sequence_files = sorted([f for f in glob.glob(directory + "/*.fastq")])
    try:
        os.mkdir(directory + "/Split_Fraction")
    except OSError:
        print ("Creation of the directory failed. Check if it already exists!")
    else:
        print ("Successfully created the directory")  
        
    for line in sequence_files:
        
        with open(directory + "/Split_Fraction/" + line.split("/")[-1].split(".")[0] + "_log.csv", "a+", newline='') as fp:
            wr = csv.writer(fp, dialect='excel')
            wr.writerow(primers.iloc[:, 0])
        
        #Iterate through each read in a given fastq file
        with open(line, "r") as handle:
            #for sequence in tqdm(list(SeqIO.parse(handle, "fastq"))[:100]):
            for sequence in tqdm(SeqIO.parse(handle, "fastq")):
                
                seq_string = str(sequence.seq)
                seq_id = ">" + str(sequence.id)
                
                #Iterate through each of the primers and check if there is a match 
                n_matches = [0,0,0,0,0,0,0,0]
                for primer in range(len(primers)):
                    
                    if m_mismatch == 0:
                    
                        if primers.iloc[:,1][primer] in seq_string:
                            n_matches[primer] = 1
                            
                            #Output each file into sub-files based on index matches
                            with open(directory + "/Split_Fraction/" + line.split("/")[-1].split(".")[0] + "_" + primers.iloc[:,0][primer] + ".fasta", "a+") as file_object:
                                file_object.write(seq_id)
                                file_object.write("\n")
                                file_object.write(seq_string)
                                file_object.write("\n")
                            
                            break 
                                
                        else:
                            pass
                            
                    else: 
                        
                        if len(regex.findall("(" + primers.iloc[:,1][primer] + "){s<=" + str(m_mismatch) + "}", seq_string)) >= 1:
                            #Record number of matches in sequence for each index 
                            n_matches[primer] = len(regex.findall("(" + primers.iloc[:,1][primer] + "){s<=" + str(m_mismatch) + "}", seq_string))
                            
                            #Output each file into sub-files based on index matches
                            with open(directory + "/Split_Fraction/" + line.split("/")[-1].split(".")[0] + "_" + primers.iloc[:,0][primer] + ".fasta", "a+") as file_object:
                                file_object.write(seq_id)
                                file_object.write("\n")
                                file_object.write(seq_string)
                                file_object.write("\n")
                        
                        else:
                            pass
                        
                
                #Write sequence with no matches to known primers to standalone sub-file
                if sum(n_matches) < 1:
                    with open(directory + "/Split_Fraction/" + line.split("/")[-1].split(".")[0] + "_unknown.fasta", "a+") as file_object:
                        file_object.write(seq_id)
                        file_object.write("\n")
                        file_object.write(seq_string)
                        file_object.write("\n")
                else:
                    pass
                
                with open(directory + "/Split_Fraction/" + line.split("/")[-1].split(".")[0] + "_log.csv", "a+", newline='') as fp:
                    wr = csv.writer(fp, dialect='excel')
                    wr.writerow(n_matches)
    
        print(str(line.split("/")[-1].split(".")[0]))

def read_count(rna, directory, n_mismatch=0):
    primer_sequence_files = sorted([f for f in glob.glob(directory + "/Split_Fraction/*.fasta")])
    read_count_df = pd.DataFrame(data = 0, index =[x.split("/")[-1].split(".")[0]for x in primer_sequence_files], columns = rna.iloc[:, 0])
 
    for fraction in primer_sequence_files:

        with open(fraction, "r") as handle:
            for sequence in tqdm(SeqIO.parse(handle, "fasta")):

                seq_string = str(sequence.seq)
                for rna_seq in range(len(rna.iloc[:,1])):
                    
                    if n_mismatch == 0:
                        if rna.iloc[:,1][rna_seq] in seq_string:
                            read_count_df.loc[fraction.split("/")[-1].split(".")[0], rna.iloc[:,0][rna_seq]] += 1
                        else:
                            pass
                    
                    else: 
                        if len(regex.findall("(" + rna.iloc[:,1][rna_seq] + "){s<=" + str(n_mismatch) + "}", seq_string)) >= 1:
                            read_count_df.loc[fraction.split("/")[-1].split(".")[0], rna.iloc[:,0][rna_seq]] += 1
                        else:
                            pass
                        
    read_count_df.to_csv(directory + "/read_counts.csv", index=True)
    
                    
def main(): 
    
    # Initiate the parser
    parser = argparse.ArgumentParser(description="Demultiplex Parameters")
    
    # Add long and short argument
    parser.add_argument("--primer", "-p", help="Tab separated file of primer sequences to separate based on fraction")
    parser.add_argument("--grna", "-g", help="Tab separated file of gRNA sequences to produce read count file")
    parser.add_argument("--directory", "-d", help="Directory of fastq files for demultiplexing")
    parser.add_argument("--mis_primer", "-m", help="Number of mismatches allowed for positive hit of primer sequence", default=0)
    parser.add_argument("--mis_rna", "-n", help="Number of mismatches allowed for positive hit of RNA read", default=0)
    # Read arguments from the command line
    args = parser.parse_args()
    
    primers = pd.read_csv(str(args.primer), sep='\t', header=None)
    primers['length'] = primers[primers.columns[1]].str.len()
    primers.sort_values('length', ascending=False, inplace=True)
    primers.reset_index(drop=True, inplace=True)
    
    rna = pd.read_csv(str(args.grna), sep='\t', header=None)
    directory = str(args.directory)
    m_mismatch = int(args.mis_primer)
    n_mismatch = int(args.mis_rna)
    
    #Demultiplex based on primer sequences and output files into new directory called Split_Fraction
    demultiplex(primers, directory, m_mismatch = m_mismatch)
    
    #Read counts for each split fraction 
    read_count(rna, directory, n_mismatch=n_mismatch)

if __name__ ==  '__main__':
    main()       
            
        
        
        
    
        
#python read_count.py 
#-p /Users/davidchen/Documents/GitHub/Sandor_David/Miniscreen05/Fw_staggered_primers.tsv 
#-g /Users/davidchen/Documents/GitHub/Sandor_David/Miniscreen05/gRNAmotifs5.tsv 
#-d /Users/davidchen/Documents/GitHub/Sandor_David/Miniscreen05/Sequence  
        
        
#python read_count.py -p /Users/davidchen/Documents/GitHub/Sandor_David/Miniscreen05/Fw_staggered_primers.tsv -g /Users/davidchen/Documents/GitHub/Sandor_David/Miniscreen05/gRNAmotifs5.tsv -d /Users/davidchen/Documents/GitHub/Sandor_David/Miniscreen05/Sequence        
        
