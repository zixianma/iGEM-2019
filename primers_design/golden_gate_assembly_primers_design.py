#define the function SPLITINTOFRAGS to split a long sequence into fragments for Golden Gate Assembly
#sequence (a str)= the input DNA sequence 
#import Bio so that we can load and parse fasta/genbank files
#import xlsxwriter so that we can output the fragment sequences to a excel file
from Bio import SeqIO
import random
import xlsxwriter
def splitIntoFrags(sequence):
    #calculate and output the length of the total input sequence in both bp and kb 
    seq_length = len(sequence)
    print('The length of the input sequence is %(length)i bp, or %(length_kb).3f kb'% {'length':seq_length,'length_kb':seq_length/1000})
   
    #convert the user input str max_num into int
    #set 2 as the default minimum # of fragments
    max_num = int(input('Please enter the maximum number of fragments:'))
    min_num = 2
    
    #create two lists for displaying the options for # of fragments & corresponding fragment lengths
    num_options = []
    length_options = []
    for i in range(max_num - min_num + 1):
        num_options.append(min_num + i)
        length_options.append(int(seq_length / num_options[i]))
     
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
    #iterate through each sequence in the list and check for invalid base
    for i in range(len(target_seq_list)):
        target_seq_list[i].upper()
        for j in range(len(target_seq_list[i])):
            if target_seq_list[i][j]!= 'A'and target_seq_list[i][j]!= 'T'and target_seq_list[i][j]!= 'C'and target_seq_list[i][j]!= 'G':
                print('The sequence you put contains invalid base, please check your file')
    return  target_seq_list      

#the function ENZYMECUTSITESELECTION takes in a target sequence (str), 
#validates the absence of the cute sites of a list of available restriction enzymes
#and returns the DNA sequence of the cut site of a usable enzyme user chose
def enzymeCutSiteSelection(target_seq):
    enzyme_dict = {'BsaI/BsaI-HF v2':'GGTCTC','BbsI/BbsI-HF':'GAAGAC','BsmBI':'CGTCTC','BfuAI/BspMI':'ACCTGC','BtgZI':'GCGATG','SapI':'GCTCTTC'}
    cut_sites = list(enzyme_dict.values())
    enzyme_choices = list(enzyme_dict.keys())
    usable_enzymes = []
    for i in range(len(cut_sites)):
        if target_seq.upper().find(cut_sites[i]) == -1:
            usable_enzymes.append(enzyme_choices[i])
    if len(usable_enzymes) > 0:
        print('The following are the enzymes you can use:')
        for i in range(len(usable_enzymes)):
            print(i,usable_enzymes[i])
        enzyme = usable_enzymes[int(input("Select the enzyme you want(enter number only):"))]
        print("You've selected the enzyme",enzyme,"The fragment sequences with the correct overhangs will be output to the demo file.")
    else:
        print("Sorry, we couldn't find a usable restriction enzyme for the sequence you input.:(")
    return enzyme_dict[enzyme]

#the function MAKEREVERSECOMPLEMENT is used for generating the 3'-5' cut site, 
#which is the complement for the forward 5'-3' cut site. However, because we are adding
#it to the end of a fragment of 5'-3' strand, we need to reverse it
def makeReverseComplement(sequence):
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

#the function ADDOVERHANGS takes in the parameter frag_list, which is obtained by parsing fasta
#and the parameter cut_site, which is returned from the ENZYMECUTSITESELECTION function
#and outputs the frag_list containing the fragments with the correct overhangs
def addOverhangs(frag_list,cut_site):
    result_list = []
    forward_cute_site = cut_site
    reverse_cut_site = makeReverseComplement(cut_site)
    random_codon = random.choice(['A','C','T','G'])
    for i in range(len(frag_list)):
        if i == len(frag_list)-1:
            start = frag_list[0][0:4]
        else:
            start = frag_list[i+1][0:4]
        new_frag = forward_cute_site + random_codon + frag_list[i] + start + random_codon + reverse_cut_site
        result_list.append(new_frag)
    return result_list

#create an excel file for storing the output
workbook = xlsxwriter.Workbook('demo.xlsx')
worksheet = workbook.add_worksheet()
target_seq_list = parseSequence('test3.fa','fasta')

#iterate through and split each sequence in the list, write the output fragment sequences in the excel file
for i in range(len(target_seq_list)):
    num,frag_list = splitIntoFrags(target_seq_list[i])
    frag_list = addOverhangs(frag_list,enzymeCutSiteSelection(target_seq_list[i]))
for k in range(num):
    worksheet.write(k,i,frag_list[k])
workbook.close()
