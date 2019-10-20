#This python file follows the primers design guide here 
#https://warwick.ac.uk/study/csde/gsp/eportfolio/directory/pg/lsujcw/gibsonguide/
from Bio import SeqIO
import random
import xlsxwriter

#define a function PARSESEQUENCE to read and validate sequences from fasta/genbank files
#file_name (a str) = the name of the input file (e.g.'test.fa')
#file_type (a str) = the type of the input file (e.g. 'fasta')
def parseSequence(file_name,file_type):
    #read a list of SequenceRecord Objects from the file
    seq_list = list(SeqIO.parse(file_name, file_type))

    #extract the actual sequences, turn them into strings and place them into a list called target_seq_list
    target_seq_list = []
    for i in range(len(seq_list)):
        length = int(len(seq_list[i].seq))
        seq = str(seq_list[i].seq[0:length])
        target_seq_list.append(seq)
    #iterate through each sequence in the list and check for invalid bases
    for i in range(len(target_seq_list)):
        target_seq_list[i].upper()
        for j in range(len(target_seq_list[i])):
            if target_seq_list[i][j]!= 'A'and target_seq_list[i][j]!= 'T'and target_seq_list[i][j]!= 'C'and target_seq_list[i][j]!= 'G':
                print('The sequence you input contains invalid base, please check your file')
    return  target_seq_list      


def splitIntoFrags(sequence):
    #calculate and output the length of the total input sequence in both bp and kb 
    seq_length = len(sequence)
    print('The length of the input sequence is %(length)i bp, or %(length_kb).3f kb'% {'length':seq_length,'length_kb':seq_length/1000})
   
    #convert the user input str max_num into int
    #set 2 as the default minimum # of fragments
    # max_num = int(input('Please enter the maximum number of fragments:'))
    # min_num = 2
    
    #create two lists for displaying the options for # of fragments & corresponding fragment lengths
    num_options = [2,3,4]
    length_options = [seq_length / i for i in num_options]

    #output options
    print('You have the following options: number of fragments: ',num_options,'corresponding mean length per fragment',length_options)

    #let user choose the # of fragments and set it to target_num
    #output the sequences of split fragments according to user input
    target_num = int(input('Select the number of fragments you want: '))
    target_length = length_options[num_options.index(target_num)]
    frag_list = []
    for i in range(target_num):
        if i == target_num - 1:
            frag = sequence[int(target_length * i):]
        else: 
            frag = sequence[int(target_length * i):int(target_length * (i+1))]
        frag_list.append(frag)
        print('fragment sequence',i+1,frag)
    return target_num,frag_list

#This function is a helper function for generating a 5'-3' reverse primer, or the 
#complement strand of a forward primer 
def generate_complement_strand(sequence):
    comp_Base_Pairs = {
    "A" : "T",
    "C" : "G",
    "T" : "A",
    "G" : "C"
    }
    complement = ''
    length = len(sequence)
    for i in range(length):
        complement = complement + comp_Base_Pairs[sequence[length - i - 1]]
    return complement

    
#This function is for calculating the melting temperature of a given priner 
#according the function Tm = 81.5 + 0.41 * gc_percentage - 675 / len(primer_sequence) - mismatch_percentage
def calculate_melting_temp(dna_sequence):
    gc_num = 0
    for i in range(len(dna_sequence)):
        if dna_sequence[i] == "G" or dna_sequence == "C":
            gc_num += 1
    gc_content = gc_num / len(dna_sequence) * 100 
    melting_temperature = 81.5 + 0.41 * gc_content - 675 / len(dna_sequence) 
    # print(melting_temperature)
    return melting_temperature

def generate_primers(frags_list,epsilon = 5):
    forward_primers,reverse_primers = [],[]
    for i in range(len(frags_list)):   
        former_frag = frags_list[i]
        latter_frag = frags_list[0] if i == len(frags_list)-1 else frags_list[i+1]

        best_former_temp_diff, best_latter_temp_diff = epsilon,epsilon
        best_former_overlap, best_latter_overlap = "", ""

        for pos in range(1,len(former_frag)):
            if former_frag[-1-pos] == 'G' or former_frag[-1-pos] == 'C':
                former_overlap = former_frag[-1-pos:-1]
                former_temp = calculate_melting_temp(former_overlap)
                if abs(former_temp - 60) <= best_former_temp_diff:
                    best_former_temp_diff = abs(former_temp - 60)
                    best_former_overlap = former_overlap
        
        
        for pos in range(1,len(latter_frag)):  
            if latter_frag[pos] == 'G' or latter_frag[pos] == 'C':
                latter_overlap = latter_frag[0:pos]
                latter_temp = calculate_melting_temp(latter_overlap)
                if abs(latter_temp - 60) <= best_latter_temp_diff:
                    best_latter_temp_diff = abs(latter_temp - 60)
                    best_latter_overlap = latter_overlap
 
        forward_primer = best_former_overlap[-21:-1] + best_latter_overlap if \
        len(best_former_overlap) >= 20 else best_former_overlap + best_latter_overlap

        best_latter_overlap_reverse = generate_complement_strand(best_latter_overlap)
        best_former_overlap_reverse = generate_complement_strand(best_former_overlap)
        reverse_primer = best_latter_overlap_reverse[-21:-1] + best_former_overlap_reverse if \
        len(best_latter_overlap_reverse) >= 20 else best_latter_overlap_reverse + best_former_overlap_reverse

        forward_primers.append(forward_primer)
        reverse_primers.append(reverse_primer)
    print("Our autogenerated former primers (5' to 3') are %s and reverse primers are %s" %(forward_primers,reverse_primers))
    return forward_primers,reverse_primers


#   former_lengths = [lambda a : len(a) for a in for potential_former_overlaps]
#         latter_lengths = [lambda a : len(a) for a in for potential_latter_overlaps]

#         print('For formerou have the following options: number of fragments: %s with the \
#         corresponding lengths %s and temperatures %s' % (potential_former_overlaps,former_lengths,former_temps))
#         print('You have the following options: number of fragments: %s with the \
#         corresponding lengths %s and temperatures %s' % (potential_former_overlaps,former_lengths,former_temps))

#     #let user choose the # of fragments and set it to target_num
#     #output the sequences of split fragments according to user input
#     target_num = int(input('Select the number of fragments you want: '))


if __name__ == "__main__":
    #create an excel file for storing the output
    workbook = xlsxwriter.Workbook('demo.xlsx')
    worksheet = workbook.add_worksheet()
    target_seq_list = parseSequence('test3.fa','fasta')

    #iterate through and split each sequence in the list, write the output fragment sequences in the excel file
    for i in range(len(target_seq_list)):
        num,frag_list = splitIntoFrags(target_seq_list[i])
        generate_primers(frag_list)
    

    workbook.close()