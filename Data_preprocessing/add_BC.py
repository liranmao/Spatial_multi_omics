from Bio.SeqIO.QualityIO import FastqGeneralIterator
from gzip import open as gzopen
import pandas as pd

import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i1", "--input_R1", required=True, help="input file R1")
ap.add_argument("-i2", "--input_R2", required=True, help="input file R2")
ap.add_argument("-l", "--input_list", required=True, help="input tagged list")

ap.add_argument("-o1", "--output_R1", required=True, help="output file R1")
ap.add_argument("-o2", "--output_R2", required=True, help="output file R2")
args = vars(ap.parse_args())

input_file_R1 = args["input_R1"]
input_file_R2 = args["input_R2"]

input_list=args["input_list"]

output_file_R1 = args["output_R1"]
output_file_R2 = args["output_R2"]

df_list = pd.read_table(input_list)
df_list = df_list.rename(columns={'#Annotation': 'Annotation'})

dict_list = dict(zip(df_list['Annotation'], zip(df_list['Match_result'], df_list['Barcode'])))

reads_count = 0

with gzopen(input_file_R1, "rt") as in_handle_R1, open(input_file_R2, "rt") as in_handle_R2, open(output_file_R1, "w") as out_handle_R1, open(output_file_R2, "w") as out_handle_R2:
    for (title1, seq1, qual1), (title2, seq2, qual2) in zip(FastqGeneralIterator(in_handle_R1), FastqGeneralIterator(in_handle_R2)): 
        if(title1.split()[0] != title2.split()[0]):
            print("Error! Two reads are not in the same order!")
            break

        if (dict_list[title2][0] == 'MATCHED_PERFECTLY') or (dict_list[title2][0] == 'MATCHED_UNAMBIGUOUSLY'):
            barcode = dict_list[title2][1]   # double check whether it is from read1 or read2
            new_title1 = barcode + ":" + title1
            new_title2 = barcode + ":" + title2
               
            out_handle_R1.write("@%s\n%s\n+\n%s\n" % (new_title1, seq1, qual1))
            out_handle_R2.write("@%s\n%s\n+\n%s\n" % (new_title2, seq2, qual2))

            reads_count += 1

    print("Total matched reads are: ", reads_count)