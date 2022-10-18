import Bio
from Bio import SeqIO
import timeit
import matplotlib.pyplot as plt
import copy
from Levenshtein import distance as lev
from termcolor import colored





# create dictionary in which the kmer is the key and the frequency of the kmer is its corresponding value
def count_kmer(file, k):
  
    beg = timeit.default_timer()

    kmer_dict = {}
    for data in Bio.SeqIO.parse(file, "fasta"):
        # select sequences read from the file
        seq = str(data.seq)
        for i in range(len(seq)-k+1):
            # divide the sequence into overlapping subsections of length k
            key = seq[i:i+k]
            if key not in kmer_dict:
                kmer_dict[key] = 0
            # increase value correpsonding to each kmer based on its frequency
            kmer_dict[key] +=1

    end1 = timeit.default_timer()
    time = end1 - beg
    
    print(f"Time required to create the dictionary : {(time)/60:.3f} minutes")
    print(f'Total Number of K-mers : {len(kmer_dict)}')

  ## -- Counting Kmers-- ##
    frequency = {}
    kmer_values = list(kmer_dict.values())
    kmer_values.sort()
    for v in kmer_values:
        if v in frequency.keys():
             frequency[v] +=1 
        else:
             frequency[v] = 1
    plot_dict = {}
    for i in range(1,1000):
        if i in frequency.keys():
            plot_dict[i] = frequency[i]
    keys = list(plot_dict.keys())
    values = list(plot_dict.values())
    results = (keys, values)
    end2 = timeit.default_timer()
    time2 = end2 - beg
    print(f"Time required to count the k-mers : {(time2)/60:.3f} minutes")
    return results,kmer_dict
#\


  

 ## -- Filter out the weak Kmers based on the emperical cut off frequency -- ##

def filter_kmer(dictionary, cf = 10):
    print(f'Initial number of the kmers: {len(dictionary)}')
    start_time = timeit.default_timer()
    for key,value in dict(dictionary).items():
        if value < cf:
            del dictionary[key]
    filter_time = timeit.default_timer()
    filtering_time = filter_time - start_time
    print("Filtering completed.!")
    print(f"It takes: {(filtering_time):.2f} seconds to filter the dictionary")
    print(f'Updated Number of k-mers: {len(dictionary)}')
    return dictionary


def genome_Assembly(kmers1: dict[str, int], kmers2: dict[str, int], k: int, sort: bool = True):
    # build reads from kmers that are only in one of the two strains
    print("SNP DETECTION ACTIVATED.......")
    start_time = timeit.default_timer()
    unique1 = [kmer for kmer in kmers1 if kmer not in kmers2]
    unique2 = [kmer for kmer in kmers2 if kmer not in kmers1]
    genome = [unique1, unique2]
    for g in genome:
        z = ""
        # z is None when no sequence can be merged
        while z is not None:
            z = None
            for x in g:
                for y in g:
                    if x != y:
                        # y stars with x
                        if x[:k - 1] == y[len(y) - k + 1:]:
                            z = y + x[k - 1:]
                            #print(z)
                        # x starts with y
                        elif x[len(x) - k + 1:] == y[:k - 1]:
                            z = x + y[k - 1:]
                            #print(z)

                        # if the two sequences are concatenated, they can be removed from the list and the resulting
                        # sequence is inserted
                        if z is not None:
                            g.remove(x)
                            g.remove(y)
                            # insert at the beginning, may be more efficient
                            g.insert(0, z)
                            #print(z)
                            break
                if z is not None:
                    break

    if sort:
        g.sort(reverse=True, key=lambda x: len(x))
        #print(only1)
        g.sort(reverse=True, key=lambda x: len(x))
        #print(only2)

    return unique1, unique2,start_time


def print_snps(t: float,sequences1: list[str], sequences2: list[str], threshold: int = 10):
    count = 0
    print()
    print("Detected SNPs:")
    for x in sequences1:
        best_dist = float('inf')
        for y in sequences2:
            dist = lev(x, y)
            if dist < best_dist:
                best_string = y
                best_dist = dist

        if best_dist < threshold:
            count += 1
            print(f"SNP {count}:")
            snp_pos = []
            for i, (c1, c2) in enumerate(zip(x, best_string)):
                if c1 != c2:
                    snp_pos.append(i)

            print("Wild")
            for i, c in enumerate(x):
                if i in snp_pos:
                    print(colored(c, 'red'), end='')
                else:
                    print(c, end='')
            print("\n---")
            print("Mutated")
            for i, c in enumerate(best_string):
                if i in snp_pos:
                    print(colored(c, 'red'), end='')
                else:
                    print(c, end='')
            print()
    end = timeit.default_timer()
    time = end - t
    print(f"It takes: {(time):.2f} seconds to detect the SNPs from the filtered K-mers")




