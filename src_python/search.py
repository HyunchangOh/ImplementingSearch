# Implementing Search in Python

#1-1 Group Name

#1-2 Naive Search Algorithm
import iv2py as iv

def import_fasta_seq_from(filename:str):
    """
    :param filename: name of the file in the 'data' directory (just the filename, not the path)
    :returns: a list of sequences of each record
    """
    records = []
    for record in iv.fasta.reader(file="../data/"+filename):
        records.append(record.seq)
    return records

def naive_search_single(substring:str, string:str):
    """
    :param substring: the substring to be searched for
    :param string: the reference string in which the substring will be searched
    :returns: list of int, that contains the starting index of all occurrences of the substring in the reference string
    """
    store = []
    while True:
        i = string.find(substring)
        if i==-1:
            break
        store.append(i)
        string=string.replace(substring,"%"*len(substring),1)
    return store

@profile
def naive_search(reference_filename:str="hg38_partial.fasta.gz", reads_filename:str="illumina_reads_100.fasta.gz", queries:int=None):
    """
    :param reference_filename: name of the file in the 'data' directory that stores the reference sequence. This should have only one record in the file.
    :param reads_filename: name of the file in the 'data' directory that stores the reads sequences. This may have very many records.
    :param queries: if entered, only the first n reads will be searched.
    :returns: list of list of int containing all the starting indices of each occurrence of the read in the reference. 
    """
    reference_sequence = import_fasta_seq_from(reference_filename)
    assert len(reference_sequence)==1, "Reference file contains more than one record"
    reference_sequence = reference_sequence[0]

    read_sequences = import_fasta_seq_from(reads_filename)

    searches = []
    if queries:
        for i in range(queries):
            read = read_sequences[i]
            searches.append(naive_search_single(read,reference_sequence))
    else:
        for read in read_sequences:
            searches.append(naive_search_single(read,reference_sequence))
    return searches

#1-3 Suffix Array Based Search
def binary_search(substring:str, suffix_array:list[int], reference:str):
    """
    :param substring: substring to be searched in the reference
    :param suffix_array: suffix array of the reference string
    :param reference: reference string in which the substring will be searched
    :returns: the index in the suffix array which contains the index where the substring occurs in the reference
    """
    right = len(suffix_array)
    left = 0
    while left<right:
        middle = (right+left)//2
        if reference[suffix_array[middle]:] < substring:
            left = middle + 1
        else:
            right = middle
    return left

import copy
def suffix_search_single(substring, suffix_array, reference):
    """
    :param substring: substring to be searched in the reference
    :param suffix_array: suffix array of the reference string
    :param reference: reference string in which the substring will be searched
    :returns: all the indices where the substring occurs in the reference
    """
    first = binary_search(substring,suffix_array,reference)

    left = first
    while left > 0:
        i = suffix_array[left]
        if reference[i: i + len(substring)] != substring:
            left+=1
            break
        else:
            left -=1

    right = first
    while right < len(suffix_array):
        i = suffix_array[right]
        if reference[i: i + len(substring)] != substring:
            break
        else:
            right +=1
    answer = copy.deepcopy(suffix_array[left:right])
    answer.sort()
    return answer

@profile
def suffix_search(reference_filename:str="hg38_partial.fasta.gz", reads_filename:str="illumina_reads_100.fasta.gz", queries:int=None):
    """
    :param reference_filename: name of the file in the 'data' directory that stores the reference sequence. This should have only one record in the file.
    :param reads_filename: name of the file in the 'data' directory that stores the reads sequences. This may have very many records.
    :param queries: if entered, only the first n reads will be searched.
    :returns: list of int containing the starting indices of each read in the reference. -1 if not found.
    """
    reference_sequence = import_fasta_seq_from(reference_filename)
    assert len(reference_sequence)==1, "Reference file contains more than one record"
    reference_sequence = reference_sequence[0]
    
    start = time.time()
    suffix_array = iv.create_suffixarray(reference_sequence)
    end = time.time()
    print("Suffix Array Constructed Time:"+str(end-start))
    read_sequences = import_fasta_seq_from(reads_filename)

    searches = []
    if queries:
        for i in range(queries):
            read = read_sequences[i]    
            searches.append(suffix_search_single(read,suffix_array,reference_sequence))
    else:
        for read in read_sequences:
            searches.append(suffix_search_single(read,suffix_array,reference_sequence))
    return searches


#1-4 Benchmark (runtime and memory) your solutions for 1’000, 10’000, 100’000 1’000’000 queries of length 100.
import time
def measure_time(query_size:int,length_size:int,suffix:bool=True, fm:bool = False):
    """
    :param query_size: number of searches made.
    :param length_size: length of the read
    :param suffix: True = Suffix array, False = naive
    :param fm: True = perform fmindex search, False = resort to suffix
    :returns: amount of time elapsed during the execution of search
    """
    start = time.time()
    if fm:
        fmindex_search(reads_filename="illumina_reads_%d.fasta.gz"%(length_size),queries=query_size)
    elif suffix:
        suffix_search(reads_filename="illumina_reads_%d.fasta.gz"%(length_size),queries=query_size)
    else:
        naive_search(reads_filename="illumina_reads_%d.fasta.gz"%(length_size),queries=query_size)
    end=time.time()
    return end-start

#1-5 Benchmark (runtime) queries of the length 40, 60, 80, and 100 with a suitable number of queries.
#Same function from 1-4 will be used for this task


#2-2 fmindex based search
@profile
def fmindex_search(reference_filename:str="hg38_partial.fasta.gz", reads_filename:str="illumina_reads_100.fasta.gz", queries:int=None):
    """
    :param reference_filename: name of the file in the 'data' directory that stores the reference sequence. This should have only one record in the file.
    :param reads_filename: name of the file in the 'data' directory that stores the reads sequences. This may have very many records.
    :param queries: if entered, only the first n reads will be searched.
    :returns: list of int containing the starting indices of each read in the reference. -1 if not found.
    """
    reference_sequence = import_fasta_seq_from(reference_filename)
    assert len(reference_sequence)==1, "Reference file contains more than one record"
    reference_sequence = reference_sequence[0]
    
    start = time.time()
    fm_index = iv.fmindex([reference_sequence,"ACTG"])
    end = time.time()
    print("FM index Constructed in Time:"+str(end-start))
    read_sequences = import_fasta_seq_from(reads_filename)

    searches = []
    if queries:
        for i in range(queries):
            read = read_sequences[i]    
            res = fm_index.search(read)
            if len(res)>0:
                for i in res:
                    searches.append([i[1]])
            else:
                searches.append([])
    else:
        for read in read_sequences:
            res = fm_index.search(read)
            if len(res)>0:
                for i in res:
                    searches.append([i[1]])
            else:
                searches.append([])
    return searches

# %load_ext memory_profiler
a = input("Welcome to ImplementSearch!!! Press Enter to Continue.")
b = input("""
    Copyright(c) Group 1: Antonio Alfaro de Prado, Eleanor Alspaugh, Hyunchang Oh
    Team "SchnappiDasKleineKrokodil"
    Choose what you want to do.
    a. Perform Task 1-4 (automatic)
    b. Perform Task 1-5 (automatic)
    c. Try Small Sized Queries (automatic, and implemented!)
""")
if b.lower()=="a":
    queries = [10000] #, 100, 1000,10000,100000]
    for q in queries:
        print("Query Size is: "+str(q))
        print("Suffix Array Search Time:")
        print(measure_time(q,100,True))
        print("-------------------")

        print("Naive Search Time:")
        print(measure_time(q,100,False))
        print("-------------------")

        print("FM Index Search Time:")
        print(measure_time(q,100,False,True))
        print("-------------------")
        print("")
if b.lower()=="b":
    lengths = [40,60,80,100]
    for l in lengths:
        print("Length Size is: "+str(l))
        print("Suffix Array Search Time:")
        print(measure_time(1000,l,True))
        print("Naive Search Time:")
        print(measure_time(1000,l,False))
        print("FM Index Search Time:")
        print(measure_time(1000,l,False,True))
if b.lower()=="c":
    print("Suffix Array Search Results: ")
    print(suffix_search(queries=20))
    print("Naive Search Results: ")
    print(naive_search(queries=20))
    print("FM index Search Results: ")
    print(fmindex_search(queries=20))