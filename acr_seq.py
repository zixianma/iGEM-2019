import numpy as np
import csv
import random

def main():
    filename = "core_dataset.csv"
    aa_seq_list = []
    data = loadCsvData(filename)
    for acr in data:
        #Can add contional statements to filter the existing protein database
        acr_type = acr[2]
        if acr_type == "AcrF7":
            aa_seq = acr[6]
            if aa_seq != "Seq":
                aa_seq_list.append(aa_seq)
    
    #Get a len_dict storing the distribution of the a.a. lengths and get the max_length 
    len_dict, max_len = getSeqLengthsDistAndMax(aa_seq_list)
    print("There are a total of %d different amino acid sequences of this protein." % (len(aa_seq_list)))
    # print("The lengths distribution of the proteins is %s." % (len_dict))
    print("The maximum length of these aa sequences is %d." % (max_len))

    #Get the dictionary aa_dist which stores the distribution of a.a.s at each position
    aa_dist = getAADistrubution(aa_seq_list)

    #Convert the counts of a.a.s into probabilities
    aa_dist = convertCountIntoProb(aa_dist)
    new_data = rearrangeDataForHeatmap(aa_dist,"acrVIB"position- 34)
    #Find out at which positions there are more than one possible a.a.s, and those are 
    #the positions of interest
    pos_of_interest =[]
    for pos in aa_dist:
        if len(aa_dist[pos]) > 1:
            pos_of_interest.append(pos)
    # print(len(pos_of_interest),pos_of_interest)
    # print(aa_dist)

    #Randomly generates a length based on the a.a. lenghts distribution 
    target_len = generateLength(len_dict)
    print("The randomly generated target length is %d." % (target_len))
    #Generate a novel a.a. sequence based on the distributions of different a.a.s
    # at different locations, assuming that they are INDEPENDENT of each other
    generated_aa_seq,length = generateAA(aa_dist,pos_of_interest,target_len)
    while generated_aa_seq in aa_seq_list:
        generated_aa_seq,length = generateAA(aa_dist,pos_of_interest,target_len)
    print("Assume that the amino acids at different positions are independently distributed, we have: %s, length: %d" % (generated_aa_seq,length))

    #Generate a novel a.a. sequence based on the distributions of different a.a.s
    # at different locations, assuming that they are DEPENDENT of each other
    generated_aa_seq_de = generateAAAsIfDependent(aa_seq_list,target_len)
    while generated_aa_seq in aa_seq_list:
        generated_aa_seq_de  = generateAAAsIfDependent(aa_seq_list,target_len)
    print("If the amino acids at different positions are dependent, we have: %s, length: %d" % (generated_aa_seq_de,len(generated_aa_seq)))

#loadCsvData() is a function for loading CSV data
def loadCsvData(fileName):
    matrix = []
    # open a file
    with open(fileName) as f:
        reader = csv.reader(f)

        # loop over each row in the file
        for row in reader:

            # cast each value to a float
            doubleRow = []
            for value in row:
                doubleRow.append(value)

            # store the row into our matrix
            matrix.append(doubleRow)
    return matrix

#Given a list of a.a. sequnces, calculate the lengths distribution and their maximum
def getSeqLengthsDistAndMax(seq_list):
    len_dict = {}
    max_len = 0
    total = len(seq_list)
    for seq in seq_list:
        #Set the maximum length
        if len(seq) > max_len:
            max_len = len(seq)
        #Add to the probabiliy of a length 
        if len(seq) in len_dict:
            len_dict[len(seq)] += 1 / total
        else: 
            #Set the key if there isn't already one
            len_dict.setdefault(len(seq)) 
            len_dict[len(seq)] = 0
            len_dict[len(seq)] += 1 / total 
    return len_dict,max_len

#This function generates a target length of the a.a. sequence we want based
#on the distribution of the a.a. lengths
def generateLength(len_dict):
    total_aa = 0
    prob = random.random()
    p_len = 0.0
    #For each possible length, get its probability and return a length according 
    #to the randomly generated number
    for length in list(len_dict.keys()):
        p_len += len_dict[length]
        if prob <= p_len:
            total_aa = length
            break
    return total_aa

#Obtain the ditribution of a.a.s at different locations given a list of a.a. sequnces
def getAADistrubution(seq_list):
    aa_dict = {}
    for seq in seq_list:
        for i in range(len(seq)):
            #Loop over each position of each sequence in the list
            if i in aa_dict:
                if seq[i] in aa_dict[i]:
                    aa_dict[i][seq[i]]+= 1 
                else:
                    #Set the key of an a.a. if there isn't already one
                    aa_dict[i].setdefault(seq[i]) 
                    aa_dict[i][seq[i]] = 0
                    aa_dict[i][seq[i]] += 1 
            else:
                #Set the key of a location if first time receiving this loc
                aa_dict.setdefault(i) 
                aa_dict[i] = {}
                #Set the key of an a.a. at a particular loc and update the count
                aa_dict[i].setdefault(seq[i]) 
                aa_dict[i][seq[i]] = 0
                aa_dict[i][seq[i]] += 1   
    return aa_dict

#This fucntion converts the counts distribution of a.a.s at different locs
#into probability distribution 
def convertCountIntoProb(aa_dict):
    for pos in aa_dict:
        total = 0
        for aa in aa_dict[pos]:
            total += aa_dict[pos][aa]
        for aa in aa_dict[pos]:
            aa_dict[pos][aa] = aa_dict[pos][aa] / total 
    return aa_dict


#This function generates a novel a.a. sequnce assuming that the a.a. are independently distributed
#across positions
def generateAA(aa_dict,pos_of_interest,total_aa):
    result = ""
    for i in range(total_aa):
        aa_list = list(aa_dict[i].keys())
        #If the position is not in our interested ones, meaning that there's only one a.a. at this 
        #position, then we just add this a.a. to our sequence, because its probability is 1
        if i not in pos_of_interest:
            aa =  aa_list[0]
            result += aa
        else:
            p = random.random()
            p_aa = 0.0
            #If we do need to choose one a.a. from multiple ones, then we choose one 
            #based on the probability distribution
            for aa in aa_list:
                p_aa += aa_dict[i][aa]
                if p <= p_aa:
                    result += aa
                    break
    return result,len(result)

#This is a helper function for generateAAAsIfDependent, and it's getting the 
#distribution of a.a.s at the next position given start index and a.a.
def getNextDist(seq_list,index,index_aa):
    #Filter the list of sequences by including only those with the specific a.a. at the 
    #specific index
    if index > 0:
        seq_list = [seq for seq in seq_list if len(seq) > index and seq[index-1] == index_aa]
    next_dist = {}
    for seq in seq_list:
        #Make sure that the starting index doesn't exceed the length of the sequence
        if seq[index] in next_dist:
            next_dist[seq[index]] += 1 
        else: 
            next_dist.setdefault(seq[index]) 
            next_dist[seq[index]] = 0
            next_dist[seq[index]] += 1 
    #Convert the counts distribution to a probability distribution
    total = 0
    for aa in next_dist:
        total += next_dist[aa]
    for aa in next_dist:
        next_dist[aa] = next_dist[aa] / total 
    # print(index,start_aa,next_dist)
    return seq_list,next_dist

#This function is for generating a.a. sequences assuming that the a.a.s are dependently
#distributed, which is most likely the case
def generateAAAsIfDependent(seq_list,aa_length):
    result = ""
    for i in range(aa_length):
        #Make sure that there's either something in "result" or there isn't a starting a.a.
        prev_aa = result[i-1] if i > 0 else ""
        #Get the distribution 
        seq_list,next_dist = getNextDist(seq_list,i,prev_aa)
        #Use the same scheme to add a.a. to the sequence
        aa_list = list(next_dist.keys())
        p = random.random()
        p_aa = 0.0
        for aa in aa_list:
            p_aa += next_dist[aa]
            if p <= p_aa:
                result += aa
                break
    return result

def rearrangeDataForHeatmap(aa_dict,acr_type):
    filename = "heatmap_data_" + acr_type + ".csv"
    writer = csv.writer(open(filename,'w'))
    writer.writerow(["group","variable","value"])
    table = []
    aa_list = []
    pos_list = list(aa_dict.keys())
    for pos in pos_list:
        for aa in list(aa_dict[pos].keys()):
            if aa not in aa_list:
                aa_list.append(aa)
    for pos in pos_list:
        for aa in aa_list:
            row = []
            row.append(pos)
            if aa in list(aa_dict[pos].keys()):
                row.append(aa)
                row.append(aa_dict[pos][aa])
            else:
                row.append(aa)
                row.append(0)
            table.append(row)
            writer.writerow(row)
    return table



# This if statement passes if this
# was the file that was executed
if __name__ == '__main__':
	main()


