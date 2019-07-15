#This is a class for desiging primers for QuickChange PCR with (1-3 nucleotides) mutations 
class PrimersDesign:
    #Initialize the class by defining the amino acid map, where one-letter 
    #a.a. code is mapped to all its corresponding codons
    def __init__(self):
        self.aa_map = {
        "F" : ["UUU","UUC"],
        "L" : ["UUA","UUG","CUU","CUC","CUA","CUG"],
        "I" : ["AUU","AUC","AUA"],
        "V" : ["GUU","GUC","GUA","GUG"],
        "S" : ["UCU","UCC","UCA","UCG","AGU","AGC"],
        "P" : ["CCU","CCC","CCA","CCG"],
        "T" : ["ACU","ACC","ACA","ACG"],
        "A" : ["GCU","GCC","GCA","GCG"],
        "Y" : ["UAU","UAC"],
        "H" : ["CAU","CAC"],
        "Q" : ["CAA","CAG"],
        "N" : ["AAU","AAC"],
        "K" : ["AAA","AAG"],
        "D" : ["GAU","GAC"],
        "E" : ["GAA","GAG"],
        "C" : ["UGU","UGC"],
        "R" : ["CGU","CGC","CGA","CGG","AGA","AGG"],
        "G" : ["GGU","GGC","GGA","GGG"],
        "W" : ["UGG"],
        "M" : ["AUG"]
        }

    #This function is a helper function for generating a 5'-3' reverse primer, or the 
    #complement strand of a forward primer 
    def generate_complement_strand(self,sequence):
        comp_Base_Pairs = {
        "A" : "U",
        "C" : "G",
        "T" : "A",
        "G" : "C",
        "U" : "A"
        }
        complement = ''
        length = len(sequence)
        for i in range(length):
            complement = complement + comp_Base_Pairs[sequence[length - i - 1]]
        return complement

    #This function is for calculating the gc content (the percentage of the bases G and C
    #of a given primer 
    def calculate_gc_content(self,dna_sequence):
        gc_num = 0
        for i in range(len(dna_sequence)):
            if dna_sequence[i] == "G" or dna_sequence == "C":
                gc_num += 1
        return gc_num / len(dna_sequence) * 100 

    #This function is for calculating the mismatch percentage of a given primer 
    def calculate_mismatch_percentage(self,original_sequence,primer_sequence):
        count = 0
        for index in range(len(original_sequence)):
            if original_sequence[index] != primer_sequence[index]:
                count += 1
        return count / len(original_sequence) * 100

    #This function is for calculating the melting temperature of a given priner 
    #according the function Tm = 81.5 + 0.41 * gc_percentage - 675 / len(primer_sequence) - mismatch_percentage
    def calculate_melting_temp(self,gc_content,primer_sequence,mismatch_percentage):
        gc_content = int(gc_content)
        melting_temperature = 81.5 + 0.41 * gc_content - 675 / len(primer_sequence) - mismatch_percentage
        return melting_temperature

    #This is the function for generating both forward and reverse primers 
    #based on the user-specified mutations and conditions 
    def generate_primers(self,pos_num,target_residue,dna_sequence,epsilon):
        #Decide the position of the mutation, and the furthermost start and end positions 
        #which are 21 nt away because the maximum length of an optimal primer is 45 nt
        mutate_pos = (pos_num - 1) * 3
        start = mutate_pos - 21
        end = mutate_pos + 21 + 3

        #Get all the possible condons of the target residue 
        target_codons = self.aa_map[target_residue]
        for_primers, rev_primers = [],[]

        #Loop over the primers of all potential lengths (25-45 nt, with the mutated 
        #codon in the middle) and only keep those that satisfy the conditions 
        for i in range(10):
            start_pos = start + i
            end_pos = end - i
            # print(start_pos,end_pos)
            #Get the part of the primer from the original sequence 
            original_sequence = dna_sequence[start_pos:end_pos]
            for j in range(len(target_codons)):
                #Modify the sequnce by swapping the middle three nts
                for_primer = dna_sequence[start_pos:mutate_pos] + target_codons[j] + dna_sequence[mutate_pos + 3:end_pos]
                gc_content = self.calculate_gc_content(for_primer)
                #print(original_sequence,for_primer,gc_content)

                #Select only those primers with a gc content greater than 40%
                if gc_content >= 30.0:
                    mismatch_percentage = self.calculate_mismatch_percentage(original_sequence,for_primer)
                    Tm = self.calculate_melting_temp(gc_content,for_primer,mismatch_percentage)
                    #print(mismatch_percentage,Tm)

                    #Select only those primers with a melting temperature around 78 degrees celcius 
                    #A small error tolerance, denoted epsilon here, can be defined by the users
                    if abs(Tm - 78) < epsilon:
                        for_primers.append(for_primer)
        #Get the reverse primers by generating the complementary strands of all the usable forward primers 
        rev_primers = [self.generate_complement_strand(for_primer) for for_primer in for_primers]
        return for_primers, rev_primers

def main():
    dna_sequence = input("Please input your DNA sequence:").upper()
    print("Got it! The length of your DNA sequence is %d" % (len(dna_sequence)))
    pos_num = int(input("At which position (as an integer between 1 and %d) do you want to mutate?"%(int(len(dna_sequence) / 3))))
    target_residue = input("What is your target protein residue (as one-letter code):").upper()
    design = PrimersDesign()
    result = design.generate_primers(pos_num,target_residue,dna_sequence,epsilon = 10)
    if len(result[0]) == 0: 
        print("You can't design any primers with the existing conditons.")
    else:
        print("Forward primers: %s; Reverse_primers%s." %(result))

if __name__ == '__main__':
    main()