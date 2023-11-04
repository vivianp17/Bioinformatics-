import sys
from Bio import SeqIO

#loading kmer files into dictionary
def loading_kmers(kmer_file):
    kmer_counts = {}
    with open(kmer_file, 'r') as kfile:
        kfile =SeqIO.parse(kmer_file, "fasta")
        for record in kfile:
            kmer_counts[str(record.seq)] = int(record.id) 
    return kmer_counts

#Reading over data file
def analzying_reads(read, kmer_counts, kmer_size, threshold):
    corrected_read = list(read)
    length = len(read)
    
    for x in range(length - kmer_size + 1):
        kmer = read[x:x+kmer_size] #slider method of kmer
        
        #If kmer meets threshold and no errors, continue loop
        if kmer in kmer_counts and kmer_counts[kmer] > threshold:
            continue
        
        #Fixing Sequence if there is Error
        else:
            for y in range(kmer_size):  
                for base in ['A', 'T', 'C', 'G']: 
                    temp_kmer = kmer[:y] + base + kmer[y+1:] #reads each nucleotide and replaces one base at a time.
                    if temp_kmer in kmer_counts and kmer_counts[temp_kmer] > threshold:                
                        corrected_read[x+y] = base
                        break                
    return ''.join(corrected_read) #joining the rest of the sequence together and returning it to main

def main():
    #variables to be taken from the command line
    kmer_count_file, data_file, kmer_size, threshold = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])

    #making the kmer dictionary
    kmer_counts = loading_kmers(kmer_count_file)

    #opening data file and outputing correctd sequences in the outfile
    with open(data_file, 'r') as data, open('outfile.fasta', 'w') as out:
        data = SeqIO.parse(data_file, "fastq")
        for record in data:
            sequence = str(record.seq)  
            corrected = analzying_reads(sequence, kmer_counts, kmer_size, threshold) #analyzing new data
            out.write('>\n' + corrected + '\n') #writing to outfile

if __name__ == '__main__':
    main()
