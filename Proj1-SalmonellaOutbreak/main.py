import argparse
import os
import time
from typing import Dict, List
import utils
import matplotlib.pyplot as plt
import numpy as np

import Bio
from Bio import SeqIO
from colorama import init
import timeit




def parse_args():
    parser = argparse.ArgumentParser(description="SNP detector")
    parser.add_argument('-f1', '--variant', type=str, help='Path to FASTA file of the variant type')

    parser.add_argument('-f2', '--wild',type=str, help='Path to FASTA file of the wild Type')      

    parser.add_argument('-k', type=int, default=31, help='Length of k-mers')

    parser.add_argument('-cf', '--cut_off_frequency', type=int,
                        help='Threshold for k-mers filtering')

    parser.add_argument('-d', '--distance-threshold', type=int, nargs='?', default=15,
                        help='Threshold for Levenshtein distance')

    parser.add_argument('-v', '--visualize', default=False, action='store_true',
                        help='Plot intermediate results')
    return parser.parse_args()




def main(args):
    k = args.k
    ## ------- 1. READ FASTA FILES -------- ##
    wild  = args.wild 
    variant = args.variant
    wild = [wild]
    variant = [variant]

    ## ------- 2. CREATE DICTIONARY FOR EACH TYPE AND COUNT THE K-MERS -------- ## 
    for f in wild :  
        print("Processing the wild type......")
        results_wild,dict_wild= utils.count_kmer(f, k)
        print("Finished Counting the k-mers for the wild type! ")

    for f in variant :  
        print("Processing the Variant type......")
        results_variant,dict_variant= utils.count_kmer(f, k)
        print("Finished Counting the k-mers for the variant type! ")

    ## ------- 3. PLOT THE K-MERS FREQUENCY HISTOGRAM -------- ## 

    if args.visualize is False :
         x,y = results_wild
         fig, ax = plt.subplots()
         
         ax.plot(x, y, label = 'k = 31')
         ax.set_yscale('log')
         ax.legend(loc = 'upper right')
         plt.title("Frequency Histogram")
         plt.ylabel('Multiplicity of k-mers')
         plt.xlabel('Frequency')
         #plt.show()
         plt.savefig('K-mer_Frequency_Histogram.pdf')


    ## ------- 4. FILTER OUT THE WEAK K-MERS -------- ## 

    if args.cut_off_frequency is None:
        thres = int(input("Select filtering threshold for wild strain: "))
    else:
         thres = args.cut_off_frequency

    print("Beginning the process of filtering out the weak k-mers in the WILD type....")
    wild_filtered = utils.filter_kmer(dict_wild, thres)
    

    print("Beginning the process of filtering out the weak k-mers in the VARIANT type....")
    variant_filtered = utils.filter_kmer(dict_variant, thres)


    ## ------- 5. MERGE K-MERS -------- ## 
    unique_wild, unique_variant,t = utils.genome_Assembly(wild_filtered, variant_filtered, k)


    ## ------- 6. OUTPUT SNPS BASED ON THE LEV DISTANCE -------- ## 
    utils.print_snps(t,unique_wild, unique_variant, args.distance_threshold)


 
  

 





    









   
            









if __name__ == '__main__':
    init()
    args = parse_args()
    main(args)
