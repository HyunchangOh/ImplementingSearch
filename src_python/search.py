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

#@profile
def naive_search(reference_filename:str="hg38_partial.fasta.gz", reads_filename:str="illumina_reads_100.fasta.gz", queries:int=None):
    """
    :param reference_filename: name of the file in the 'data' directory that stores the reference sequence. This can have more than one record in the file.
    :param reads_filename: name of the file in the 'data' directory that stores the reads sequences. This may have very many records.
    :param queries: if entered, only the first n reads will be searched.
    :returns: list of list of int containing all the starting indices of each occurrence of the read in the reference. 
    """
    # Importing reference sequence, and if it's got more than one record, joining them all together, as in the case for the whole human genome
    reference_sequence = import_fasta_seq_from(reference_filename)
    if (len(reference_sequence) != 1):
        print("Reference file contains more than one record")
        reference_sequence = ''.join(reference_sequence)
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
    :param reference_filename: name of the file in the 'data' directory that stores the reference sequence. This can have more than one record in the file.
    :param reads_filename: name of the file in the 'data' directory that stores the reads sequences. This may have very many records.
    :param queries: if entered, only the first n reads will be searched.
    :returns: list of int containing the starting indices of each read in the reference. -1 if not found.
    """
    # Importing reference sequence, and if it's got more than one record, joining them all together, as in the case for the whole human genome
    reference_sequence = import_fasta_seq_from(reference_filename)
    if (len(reference_sequence) != 1):
        print("Reference file contains more than one record")
        reference_sequence = ''.join(reference_sequence)
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
def measure_time(query_size:int, length_size:int, search:str="s"):
    """
    :param query_size: number of searches made.
    :param length_size: length of the read
    :param search: "s" = Suffix array, "n" = Naive search, "f" = FM index based search, "p" = Filtered PEX algorithm based search
    :returns: amount of time elapsed during the execution of search
    """
    start = time.time()
    if search == "s":
        suffix_search(reads_filename="illumina_reads_%d.fasta.gz"%(length_size),queries=query_size)
    if search == "n":
        naive_search(reads_filename="illumina_reads_%d.fasta.gz"%(length_size),queries=query_size)
    if search == "f":
        fmindex_search(reads_filename="illumina_reads_%d.fasta.gz"%(length_size),queries=query_size)
    if search == "p":
        filt_search(reads_filename="illumina_reads_%d.fasta.gz"%(length_size),queries=query_size)
    else:
        print("No search algorithm matches input.")
    end=time.time()
    return end-start

#1-5 Benchmark (runtime) queries of the length 40, 60, 80, and 100 with a suitable number of queries.
#Same function from 1-4 will be used for this task


#2-2 fmindex based search
@profile
def fmindex_search(reference_filename:str="hg38_partial.fasta.gz", reads_filename:str="illumina_reads_100.fasta.gz", queries:int=None):
    """
    :param reference_filename: name of the file in the 'data' directory that stores the reference sequence. This can have more than one record in the file.
    :param reads_filename: name of the file in the 'data' directory that stores the reads sequences. This may have very many records.
    :param queries: if entered, only the first n reads will be searched.
    :returns: list of int containing the starting indices of each read in the reference. -1 if not found.
    """
    # Importing reference sequence, and if it's got more than one record, joining them all together, as in the case for the whole human genome
    reference_sequence = import_fasta_seq_from(reference_filename)
    if (len(reference_sequence) != 1):
        print("Reference file contains more than one record")
        reference_sequence = ''.join(reference_sequence)
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

# Function to split a string into error+1 number of parts
def splitting(input_string:str, num_parts:int):
    """
    :param input_string: string which we want to split
    :param num_parts: number of parts in which we wish to split the string
    :returns: list which contains the parts (strings) of the original string
    """
    # Calculate the length of each part
    part_length = len(input_string) // num_parts
    remainder = len(input_string) % num_parts
    # Iterate over the number of parts
    parts = []
    start_index = 0
    for i in range(num_parts):
        # Calculate the end index for the current part
        end_index = start_index + part_length + (1 if i < remainder else 0)
        # Append the currrent part to the list
        parts.append(input_string[start_index:end_index])
        # Update the start index for the next iteration
        start_index = end_index
    return parts

# Function to verify hits, based on the number of errors allowed, from a list of lists of positions (1st dim = no of parts, 2nd dim = locations)
def verify(text:str, parts:list, parts_loc:list, error:int):
    """
    :param text: original reference text
    :param parts: parts of the original substring we wanted to search for in the original text
    :param parts_loc: list of i parts, in which in each i we find a list of locations inside the original text
    :param error: maximum number of mismatches allowed
    :param num_parts: number of parts in which we wish to split the string
    :returns: list which contains the verified locations, in order and no duplicates, of the original substring made up of the parts of var parts
    """
    loc_list = []
    # For every part
    for i in range(len(parts)):
        # Number of characters of the string part, to the left of it and to the right
        le = len(parts[i])
        L_char = len(''.join(parts[:i]))
        R_char = len(''.join(parts[(i+1):]))
        # We check the other parts
        for j in range(len(parts_loc[i])):
            # Here we get the reference position, and the left and right text to the reference position
            ref_pos = parts_loc[i][j]
            L_text = text[(ref_pos-L_char):ref_pos]
            R_text = text[(ref_pos+le):(ref_pos+R_char)]
            # Here we get the number of discrepancies to either side
            discrepanciesL = sum(a != b for a, b in zip(L_text, ''.join(parts[:i])))
            discrepanciesR = sum(a != b for a, b in zip(R_text, ''.join(parts[(i+1):])))
            if (discrepanciesL+discrepanciesR <= error):
                loc_list.append(parts_loc[i][j])
    # We return a list of positions, only those that contan discrepancies <= error, sort it, and remove duplicates
    loc_list.sort()
    loc_list = list(dict.fromkeys(loc_list))
    return loc_list

# Function to look for a substraing in a text, using filtered search, and with a maximum amount of errors allowed
def filt_search_single(text:str, s_read:str, index, nerror:int):
    """
    :param text: original reference text
    :param s_read: substring to look for in the text
    :param index: fm index to use as base for a search algorithm
    :param nerror: maximum number of mismatches allowed
    :returns: list which contains the verified locations, in order and no duplicates, of the original substring made up of the parts of var parts
    """
    read_parts = splitting(s_read, nerror + 1)
    # We look for each part of the substring individually, and store each list of locations in rparts_loc, it being a list of lists of positions
    rparts_loc = []
    for i in range(len(read_parts)):
        rparts_loc.append([])
        rparts_loc[i].append(index.search(read_parts[i]))
    # Format, we take out the second dimension which is useless, and all 0 from position lists
    for i in range(len(rparts_loc)):
        rparts_loc[i] = rparts_loc[i][0]
        for j in range(len(rparts_loc[i])):
            rparts_loc[i][j] = rparts_loc[i][j][1]
    # We retrieve the list of positions
    searches = verify(text, read_parts, rparts_loc, nerror)
    return searches

# Fast filtering algorithm
@profile
def filt_search(reference_str:str="hg38_partial.fasta.gz", reads_filename:str="illumina_reads_100.fasta.gz", queries:int=None, error:int=1):
    """
    :param reference_filename: name of the file in the 'data' directory that stores the reference sequence. This can have more than one record in the file.
    :param reads_filename: name of the file in the 'data' directory that stores the reads sequences. This may have very many records.
    :param queries: if entered, only the first n reads will be searched.
    :param error: maximum number of errors allowed in string search. 1 by default.
    :returns: list of int containing the starting indices of each read in the reference. -1 if not found.
    """
    # Importing reference sequence, and if it's got more than one record, joining them all together, as in the case for the reference human genome
    reference_sequence = import_fasta_seq_from(reference_str)
    if (len(reference_sequence) != 1):
        print("Reference file contains more than one record")
        reference_sequence = ''.join(reference_sequence)
    reference_sequence = reference_sequence[0]
    
    # Importing sequences to search for
    read_sequences = import_fasta_seq_from(reads_filename)

    # FM Index
    start1 = time.time()
    fm_index = iv.fmindex([reference_sequence,"ACTG"])
    end1 = time.time()
    print("FM index Constructed in Time:" + str(end1 - start1))

    searches = []
    if queries:
        for i in range(queries):
            read = read_sequences[i]
            res = filt_search_single(reference_sequence, read, fm_index, error)
            if len(res)>0:
                for i in res:
                    searches.append([i])
            else:
                searches.append([])
    else:
        for read in read_sequences:
            res = filt_search_single(reference_sequence, read, fm_index, error)
            if len(res)>0:
                for i in res:
                    searches.append([i])
            else:
                searches.append([])
    return searches


%load_ext memory_profiler
a = input("Welcome to ImplementSearch! Press Enter to Continue.")
b = input("""
    Copyright(c) Group 1: Antonio Alfaro de Prado, Eleanor Alspaugh, Hyunchang Oh
    Team "SchnappiDasKleineKrokodil"
    Choose what you want to do.
    a. Perform Task 1-4 (automatic)
    b. Perform Task 1-5 (automatic)
    c. Try Small Sized Queries (automatic, and implemented!)
    d. Assignment 7 Search with errors based on pigeon hole principle (demo)
""")
if b.lower()=="a":
    queries = [10000] #, 100, 1000,10000,100000]
    for q in queries:
        print("Query Size is: "+str(q))
        print("Suffix Array Search Time:")
        print(measure_time(q,100,"s"))
        print("-------------------")

        print("Naive Search Time:")
        print(measure_time(q,100,"n"))
        print("-------------------")

        print("FM Index Search Time:")
        print(measure_time(q,100,"f"))
        print("-------------------")
        print("")
if b.lower()=="b":
    lengths = [40,60,80,100]
    for l in lengths:
        print("Length Size is: "+str(l))
        print("Suffix Array Search Time:")
        print(measure_time(1000,l,"s"))
        print("Naive Search Time:")
        print(measure_time(1000,l,"n"))
        print("FM Index Search Time:")
        print(measure_time(1000,l,"f"))
if b.lower()=="c":
    print("Suffix Array Search Results: ")
    print(suffix_search(queries=20))
    print("Naive Search Results: ")
    print(naive_search(queries=20))
    print("FM index Search Results: ")
    print(fmindex_search(queries=20))
if b.lower()=="d":
    print("Demo of filtered search:")
    print(measure_time(100, 100, "p"))