import Bio 
from Bio.Blast import NCBIWWW
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
fasta_file = open(r'C:\Users\PowerEill\Downloads\dna2.fasta')
fasta_file_parse = list(SeqIO.FastaIO.FastaIterator(fasta_file))

print(f"Number of sequences in the file: {len(fasta_file_parse)}")


longest_len = max(fasta_file_parse, key = lambda x: len(x.seq))
shortest_len =min(fasta_file_parse, key = lambda x: len(x.seq))
print(f"Longest sequence length: {len(longest_len)}")
print(f"Shortest sequence length: {len(shortest_len)}")

shortest_sequence=[short for short in fasta_file_parse if len(short.seq) == len(shortest_len)]
longest_sequence=[long for long in fasta_file_parse if len(long.seq) == len(longest_len.seq)]
print(f"Number of sequences with the longest length: {len(longest_sequence)}")
print(f"Number of sequences with the shortest length: {len(shortest_sequence)}")

def find_orfs(frame,fasta_file):
    start_codon = 'ATG'   
    stop_codon = ['TAA','TAG','TGA']
    read_codon = []
    start_index = None
    for Seq in fasta_file:
        sequence = Seq.seq.upper()
        seq_len = len(sequence)

        i = frame
        while i< seq_len-2:
            codon = sequence[i:i+3]
            if codon == start_codon:
                start_index = i
                for j in range(start_index, seq_len-2, 3):
                    codon = sequence[j:j+3]
                    if codon in stop_codon:
                        stop_index = j+3
                        seq_length = stop_index - start_index
                        read_codon.append((Seq.id, start_index, stop_index, seq_length,
                                           sequence[start_index:stop_index]))
                        break
            i += 3
    return read_codon

frame1 = find_orfs(0,fasta_file_parse)
frame2 = find_orfs(1,fasta_file_parse)
frame3 = find_orfs(2,fasta_file_parse)

reverse_file = []
for seq_record in fasta_file_parse:
    rev_comp_seq = seq_record.seq.reverse_complement()
    rev_comp_record = SeqRecord(
        rev_comp_seq,  
        id=seq_record.id,  
        description="Reverse complement of " + seq_record.description)
    reverse_file.append(rev_comp_record)

comp_frame1 = find_orfs(0,reverse_file)
comp_frame2 = find_orfs(1,reverse_file)
comp_frame3 = find_orfs(2,reverse_file)

print(max(frame2,key = lambda x: x[3]))
print(max(frame3,key = lambda x: x[3]))

longest_target_ORF = [seqe for seqe in (frame1+frame2+frame3) if seqe[0] == 'gi|142022655|gb|EQ086233.1|16']
print(max(longest_target_ORF, key = lambda x: x[3]))

from collections import Counter 
count = Counter()         
for seqe in (fasta_file_parse):
    sequence = seqe.seq
    for i in range(len(sequence)-11):
        pattern = sequence[i:i+12]
        count[pattern] += 1
most_frequent_kmer, max_count = count.most_common(1)[0]
count.most_common(10)
print(f"Most frequent 6-mer: {most_frequent_kmer}, most count: {max_count}")

fasta_file_parse
