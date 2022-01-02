#!/usr/bin/env python3.8

import Bioinfo
import gzip
import matplotlib
import matplotlib.pyplot as plt
import argparse

#index file
#index_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"

# # real files
#read_file_1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
#index_file_1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
#index_file_2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
#read_file_2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"

# test files
# read_file_1 = "/home/jjacobso/bgmp/bioinformatics/Bi622/Demultiplexing-Bi622/Assignment-the-third/test_R1.fq"
# index_file_1 = "/home/jjacobso/bgmp/bioinformatics/Bi622/Demultiplexing-Bi622/Assignment-the-third/test_R2.fq"
# index_file_2 = "/home/jjacobso/bgmp/bioinformatics/Bi622/Demultiplexing-Bi622/Assignment-the-third/test_R3.fq"
# read_file_2 = "/home/jjacobso/bgmp/bioinformatics/Bi622/Demultiplexing-Bi622/Assignment-the-third/test_R4.fq"


# Argparse file option 
def read_and_index_file_input():
    parser = argparse.ArgumentParser(description = "Argparse")
    parser.add_argument ("-f", help="Read file 1 path location", required =True, type = str)
    parser.add_argument ("-s", help="Index file 1 path location", required =True, type = str)
    parser.add_argument ("-t", help="Index file 2 path location", required =True, type = str)
    parser.add_argument ("-l", help="Read file 2 path location", required =True, type = str)
    parser.add_argument ("-i", help="Index file path location", required =True, type = str)
    return parser.parse_args()

# args = read_and_index_file_input()
read_file_1 = args.f
index_file_1 = args.s
index_file_2 = args.t
read_file_2 = args.l
index_file = args.i


with open(index_file, "r") as fh:
    count = 0
    barcode_list = []
    barcode_name = []
    barcode_dictionary = {}
    barcode_counter_dictionary = {}
    for line in fh:
        line = line.split()
        if line[0] != "sample":
            barcode_list.append(line[4])
            barcode_name.append(line[3])
            barcode_dictionary[barcode_name[count]] = barcode_list[count]
            barcode_counter_dictionary[barcode_name[count]] = 0
            count +=1
        
#two copies of the barcode name list. These will be overwritten with 'file pointers' below        
barcode_name_copy1 = barcode_name.copy()
barcode_name_copy2 = barcode_name.copy()

def reverse_compliment_check(barcode_1: str, barcode_2: str):
    '''This function will take the barcodes from two files. 
   Return 1: Reverse Compliment and proper barcodes
   Return 2: Not Reverse Barcode
   Return 3: Reverse Compliment but not in Dictionary
   Return 4: Not reverse compliment but in dictionary (index hopped)
    '''
    barcode_2_RC = ""
    barcode_1_remake = ""
    reverse_direction_counter = 0
    for chr in barcode_2:
        reverse_direction_counter -= 1

        if barcode_2[reverse_direction_counter] == "A":
            barcode_2_RC = barcode_2_RC + "T"
        elif barcode_2[reverse_direction_counter]== "T":
            barcode_2_RC = barcode_2_RC + "A"
        elif barcode_2[reverse_direction_counter] == "C":
            barcode_2_RC = barcode_2_RC + "G"
        elif barcode_2[reverse_direction_counter] == "G":
            barcode_2_RC = barcode_2_RC + "C"
        elif barcode_2[reverse_direction_counter] == "N":
            barcode_2_RC = barcode_2_RC + "N"
    
    for chr in barcode_1:
        if chr == "A":
            barcode_1_remake = barcode_1_remake + "A"
        elif chr == "T":
            barcode_1_remake = barcode_1_remake + "T"
        elif chr == "C":
            barcode_1_remake = barcode_1_remake + "C"
        elif chr == "G":
            barcode_1_remake = barcode_1_remake + "G"
        elif chr == "N":
            barcode_1_remake = barcode_1_remake + "N"

    if barcode_2_RC == barcode_1_remake and barcode_1_remake in barcode_dictionary.values() and barcode_2_RC in barcode_dictionary.values():
        # 1 = reverse compliment and a proper barcode
        return 1
    if barcode_2_RC in barcode_dictionary.values() and barcode_1_remake in barcode_dictionary.values() and (barcode_2_RC != barcode_1_remake):
        # 4 = not reverse compliment but both are in dictionary
        #these are index hopped!
        return 4
    if barcode_2_RC !=  barcode_1_remake:
        # 2 = not reverse compliment
        return 2
    if barcode_2_RC == barcode_1_remake and (barcode_1_remake not in barcode_dictionary.values() and barcode_2_RC not in barcode_dictionary.values()):
        #3 = reverse compliment but not in dictionary
        return 3
    

def sort_low_quality_barcodes(Qscore_1: str, Qscore_2: str):
    '''
    Input is barcode Qscores. This converts to phred scores and 
    if the average is below Q30 it returns 6. If the average is greater or equal to 30
    it returns 5.
    '''
    Q1 = Bioinfo.qual_score(Qscore_1)
    Q2 = Bioinfo.qual_score(Qscore_2)
    if Q1 >= 30 and Q2 >= 30:
        return 5
    else: 
        return 6
        
#open all files

#create all new files R1
z =0
for barcode in barcode_name_copy1:
    pfix = barcode +"_R1"
    barcode_name_copy1[z] = open(pfix + "_output.fq", "w")
    z+=1

#create all new files R2
x =0
for barcode in barcode_name_copy2:
    pfix = barcode +"_R2"
    #print(banana)
    barcode_name_copy2[x] = open(pfix + "_output.fq", "w")
    #print(barcode_name_copy1[x])
    x+=1

#create special case files
mismatch_read1 = open("mismatch_read1.fq", "w")
mismatch_read2 = open("mismatch_read2.fq", "w")
non_existing_read1 = open("non_existing_read1.fq", "w")
non_existing_read2 = open("non_existing_read2.fq", "w")


# open real files
index_file_1 = gzip.open(index_file_1, "rt")
index_file_2 = gzip.open(index_file_2, "rt")
read_file_1 = gzip.open(read_file_1, "rt")
read_file_2 = gzip.open(read_file_2, "rt")

#open test files
# index_file_1 = open(index_file_1, "r")
# index_file_2 = open(index_file_2, "r")
# read_file_1 = open(read_file_1, "r")
# read_file_2 = open(read_file_2, "r")


#Begin fun stuff here!

#here are a bunch of counters, mostly used for debugging
my_count = 0
mod_0 = 0
mod_1 = 0
mod_2 = 0
mod_3 =0
mod_sum = 0

#here are counters for each file type
correct_barcode =0
unmatched_barcode =0
unreal_barcode = 0
index_hopped_barcode =0

#Loop through R1,R2,R3,R4
modular = 0
for I1,I2,R1,R2 in zip(index_file_1, index_file_2, read_file_1, read_file_2):
    modular +=1

    #line1
    if modular % 4 == 1:
        read1_line1 = R1
        read2_line1 = R2
        mod_1 +=1

    #line2
    if modular % 4 == 2:
        read1_line2 = R1
        read2_line2 = R2
        index1_barcode = I1
        index2_barcode = I2
        mod_2 +=1
    
    #line3
    if modular % 4 == 3:
        read1_line3 = R1
        read2_line3 = R2
        mod_3 +=1

    #line4
    if modular % 4 == 0:
        index1_qscore = I1
        index2_qscore = I2
        read1_line4 = R1
        read2_line4 = R2
        mod_0 +=1

        header1 = read1_line1 + ":" + index1_barcode
        header1 = "".join(header1.split())
        header2 = read2_line1 + ":" + index2_barcode
        header2 = "".join(header2.split())

        #Filter low quality reads out
        are_barcodes_high_quality = sort_low_quality_barcodes(index1_qscore, index2_qscore)
        if are_barcodes_high_quality == 6:
            #sort to non_existing_reads files
            
            non_existing_read1.write(header1+ "\n" + read1_line2 + read1_line3 + read1_line4)
            non_existing_read2.write(header2+ "\n" + read2_line2 + read2_line3 + read2_line4)

            unreal_barcode += 1

        #this continues to the next filter
        elif are_barcodes_high_quality == 5:
            
            #filter out the rest. 3 options
            are_barcodes_reverse_compliments = reverse_compliment_check(index1_barcode, index2_barcode)
            if are_barcodes_reverse_compliments == 2:
                #sort to non_existant files
                
                non_existing_read1.write(header1+ "\n" + read1_line2 + read1_line3 + read1_line4)
                non_existing_read2.write(header2+ "\n" + read2_line2 + read2_line3 + read2_line4)

                unreal_barcode += 1
            
            elif are_barcodes_reverse_compliments ==4:
                #index hopped and move to junk
            
                mismatch_read1.write(header1+ "\n" + read1_line2 + read1_line3 + read1_line4)
                mismatch_read2.write(header2+ "\n" + read2_line2 + read2_line3 + read2_line4)

                index_hopped_barcode +=1
                
            elif are_barcodes_reverse_compliments == 3:
                #sort to non_existing_barcodes (same file as low quality reads)
            
                non_existing_read1.write(header1+ "\n" + read1_line2 + read1_line3 + read1_line4)
                non_existing_read2.write(header2+ "\n" + read2_line2 + read2_line3 + read2_line4)

                unreal_barcode +=1

            #continue to the final filter - all these have correct barcodes
            elif are_barcodes_reverse_compliments ==1:

                 index1_barcode= "".join(index1_barcode.split())

                 for value in barcode_dictionary.keys():
                    if barcode_dictionary[value] == index1_barcode:
                        prefix = value
                 

                 prefix1 = prefix + "_R1"
                 prefix2 = prefix + "_R2"

                 y=-1
                 #here is where all the read are sent to their proper files
                 for barcode in barcode_name:
                    coconut = barcode +"_R1"
                    y +=1
                    if prefix1 == coconut:
                        barcode_name_copy1[y].write(header1 + "\n" + read1_line2 + read1_line3 + read1_line4)
                        barcode_name_copy2[y].write(header2 + "\n" + read2_line2 + read2_line3 + read2_line4)
                        barcode_counter_dictionary[barcode] +=1
                 y=-1

                 correct_barcode +=1
        #debug line         
        else:
            print("This should never appear")


#final counts exported to a file
final_counts =open("final_counts.txt", "w")
final_counts.write("Index hopped Barcodes: " + str(index_hopped_barcode) + "\n")
final_counts.write("Low quality / Nonexistant barcodes: " + str(unreal_barcode)+ "\n")
final_counts.write("Correctly matched barcodes: " + str(correct_barcode) + "\n")

#final counts + extra printed or sent to slurm file
print(" final sum mod = " , mod_1 +mod_2 + mod_3 + mod_0)
print("final count = ", correct_barcode + unreal_barcode + index_hopped_barcode)
print("Nonexistant = ", unreal_barcode)
print("Correct = ", correct_barcode)
print("Index hopped: ", index_hopped_barcode)
print("Barcode counter dictionary: ", barcode_counter_dictionary)

#Chart stuff
correct_sum = 0
for item in barcode_counter_dictionary.values():
    correct_sum += item

percent_list_barcodes = []
percent_list_barcodes_names = []

for value in barcode_counter_dictionary.values():
    percent_list_barcodes.append(value /correct_sum *100)

for key in barcode_counter_dictionary.keys():
    percent_list_barcodes_names.append(key)

#yes this could also be range(24)
num_of_barcodes = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
barcode_png_thing = list(barcode_counter_dictionary.values())

#correct barcode distribution
plt.bar(num_of_barcodes, percent_list_barcodes)
plt.title('Distribution of Correct Barcodes')
plt.xticks(num_of_barcodes, percent_list_barcodes_names, rotation = "vertical")
plt.xlabel('Barcode')
plt.ylabel('Percentage')
plt.grid(True)
plt.savefig('Individual_correct_barcodes.png')
plt.show
plt. clf()

sum_stuff = correct_barcode+index_hopped_barcode+unreal_barcode
barcode_type_counts =[(correct_barcode/sum_stuff *100), (index_hopped_barcode/sum_stuff *100), (unreal_barcode/sum_stuff *100)]
barcode_types = ["Correctly Barcoded", "Index Hopped Barcodes", "Junk Barcodes"]

#Totals barcode distribution
plt.bar([1,15,29], barcode_type_counts, width = 4)
plt.title('Distribution of Barcode Events')
plt.xticks([1,15,29], barcode_types, rotation = "horizontal")
plt.xlabel('Barcode Types')
plt.ylabel('Percentage')
plt.grid(True)
plt.savefig('Barcode_types.png')


#Close all files
mismatch_read1.close
mismatch_read2.close
non_existing_read1.close
non_existing_read2.close

index_file_1.close
index_file_2.close
read_file_1.close
read_file_2.close

#Close all files
q=0
for barcode in barcode_name_copy1:
    barcode_name_copy1[q].close
    q+=1

w =0
for barcode in barcode_name_copy2:
    barcode_name_copy2[w].close
    w+=1
