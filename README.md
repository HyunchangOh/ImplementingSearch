# ImplementingSearch
### Team "SchnappiDasKleineKrokodil"

Group 1: Antonio Alfaro de Prado, Eleanor Alspaugh, Hyunchang Oh

This is an exercise to demonstrate the power of Suffix-Array and FM-Index based searching.
Visit our Repository at https://github.com/HyunchangOh/ImplementingSearch

This implementation was developed and tested on WSL(Ubuntu) on Windows 10. 

# Python Implementation
C++ follows below.
## Dependencies
1. This program runs on Python. In order to install python, please refer to https://www.python.org/downloads/
2. This program requires installing the library iv2py. Refer to https://pypi.org/project/iv2py/ for installation.

## How to Run
In the src_python directory, simply run "python search.py" and command-line interface will pop up.

Press Enter after Welcome message.

There are three modes:
a. Perform Task 1-4
b. Perform Task 1-5
c. Perform small-scale search (20 queries) and displays the results

Choose the task you wish to perform. All results will be printed in command line.



## Implementation Description
This program has two ways of search, namely: naive search and suffix search.

Both searches make use of the 'single' version of their own (suffix_search_single() and naive_search_single()),  which searches for a single substring in the reference string. Upon calling the search function, this 'single' function will be called for each read, until all specified reads have been searched.

### Naive Search
1. The Python native function 'find' is used to search for the leftmost occurrence of the substring.
2. The leftmost occurrence will be replaced by a placeholder string.
3. The search will be called again. With the leftmost occurrence now replaced, the second occurrence will be returned by 'find', if any.
4. The search terminates if find returns -1.

### Suffix Search
1. Suffix Array will be constructed by the library function iv2py.construct_suffixarray
2. Binary search the suffix array to find the 'middle' occurrence of the substring.
3. Travel upward and downward in the suffixarray to find extra occurrences of the substring.

### Fm Index Search
1. Construct Fm Index with iv2py library.
2. Perform search and reorganize format to comply with other search implementations.

### Pigeonhole search
1. Build fm index for future searches.
2. Split substring into error+1 parts.
3. Search for parts individually with fm index based search.
4. Verify each hit individually by looking at the hit position in the text and comparing the rest of the text with the substring.
5. Save verified positions, allowing for the determined number of errors. Delete duplicate positions.

## Benchmark Comparison
![image](https://github.com/HyunchangOh/ImplementingSearch/assets/42934606/4499b395-5075-454a-a66e-ff517032d648)


For query size of 10,100,1000, the memory usage and computation time of each search method was compared, using python's native time library for computation time and memory-profiler for memory usage.
For detailed use of memory-profiler library, visit https://pypi.org/project/memory-profiler/. Memory is in MB and time is in seconds.


### Memory Usage
1. Naive search does not generate any additional data structures, so its memory use is the lowest, while Suffix Array uses significantly more memory than FM index.

### Computation Time
1. Of the three algorithms tested between assignments 1 and 2, the FMindex search appears to be the most efficient. One reason for this is its average-case linear time complexity. The FMindex search algorithm also uses a log-based time complexity for encoding of the test while suffix array is not log-based; this makes it faster than suffix array for individual searches. It is based on the Burrows Wheeler Transform (BWT) and suffix array. The FMindex search also has lower time requirements than naive search. The naive search does not require the same pre-processing that FMindex and suffix array require but this pre-processing step reduces downstream processing time during the search.
2. Pigeonhole Search time increased with more allowed maximum errors, espeically with higher query number.

# C++ Program
## Dependencies
1. This program runs with Seqan3 and libsufsort, which is already installed in 'lib' subdirectory.
2. The dependencies, especially seqan3, require specific versions of compilers (and possibly others). Check troubleshooting if you run into errors upon running 'cmake'.

## Implementation Description
### Naive Search
1. Naively iterate through the target sequence.
2. In the for-loop in 1, iterate through the substring sequence.
3. Compare character by character.

### Suffix Search
1. Build suffixarray sort with the help of libsufsort library.
2. Perform binary search on the suffix array.

### Fm Index Search
1. Construct FM index by fmindex_construct.cpp with the help of skeleton codes from SGSSGene.
2. Perform FM index search with the help of Seqan3 Library's 'search' function. (probably why this is super fast compared to other search methods..)

### Pigeonhole Search
1. Build fm index for future searches.
2. Split substring into error+1 parts.
3. Search for parts individually with fm index based search.
4. Verify each hit individually by looking at the hit position in the text and comparing the rest of the text with the substring.
5. Save verified positions, allowing for the determined number of errors. Delete duplicate positions.


## How to Build
```
$ # We are assuming you are in the terminal/console inside the repository folder
$ mkdir build # creates a folder for our build system
$ cd build
$ cmake ..    # configures our build system
$ make        # builds our software, repeat this command to recompile your software

$ ./bin/naive_search --reference ../data/hg38_partial.fasta --query ../data/illumina_reads_40.fasta --query_ct 100 # calls the code in src/naive_search.cpp

$ ./bin/suffixarray_search --reference ../data/hg38_partial.fasta --query ../data/illumina_reads_40.fasta # calls the code in src/suffixarray_search.cpp

$ ./bin/fmindex_construct --reference ../data/hg38_partial.fasta --index myIndex.index # creates an index, see src/fmindex_construct.cpp

$ ./bin/fmindex_search --index myIndex.index --query ../data/illumina_reads_40.fasta --query_ct 100 --errors 0  # searches by using the fmindex, see src/fmindex_search.cpp

$ ./bin/fmindex_pigeon_search --reference ../data/hg38_partial.fasta --index myIndex.index --query ../data/illumina_reads_40.fasta --query_ct 100 --errors 0  # searches by using the fmindex, see src/fmindex_pigeon_search.cpp
```
## Benchmark Comparison 
![image](https://github.com/HyunchangOh/ImplementingSearch/assets/42934606/768c040f-f6d7-422b-82a8-33830ba443dd)\
All running time measured in nanoseconds.
1. Running time: naive > pigeonhole > fmindex > suffix array
2. Query Substring Size: running time of naive algorithm grows linearly with substring size, logarithmically with pigeonhole, but has no impact on fmindex and suffix array search.
3. Query Size: running time of all, but fmindex, grows linearly with query size. (presumably also grows linearly with fmindex, but fmindex was super fast that this was unmeasurable.)
Note that analysis of C++ algorithm was done in a manner to explore a different aspect than its Python counterpart, but all other functions, such as changing allowed maximum error in C++ or changing query substring size in Python are all implemented and available for more exploration.

## Troubleshooting
### Pushing a Large File to a Github Repo
If pushing to a github repository returns:
fatal: the remote end hung up unexpectedly

Try (but it is not recommended to push a large file like .fasta file directly to a github repo):
```
git config http.postBuffer 524288000
```
### Seqan3 Compiler Dependency Error upon cmake
Use these commands to check/update your C, C++ compiler
```
$ gcc --version
$ g++ --version
$ which gcc #locate the gcc file your system is using
$ which g++ #locate the g++ file your system is using
$ sudo apt install gcc-12 # WARNING: seqan3 is updated regularly, and may require a different version
$ sudo apt install g++-12 # WARNING: seqan3 is updated regularly, and may require a different version
$ cd /usr/bin
$ mv gcc-12 gcc # WARNING: this will overwrite your original gcc with gcc-12 and ruin all your other programs that run with original gcc, but this is guaranteed to make seqan3 work.
$ mv g++-12 g++ # WARNING: this will overwrite your original g++ with g++-12 and ruin all your other programs that run with original gcc, but this is guaranteed to make seqan3 work.
```

### Seqan3 Unable to Locate GZIP
```
$ sudo install zlib1g
```
This installs ZLIB, but Seqan3 may still not be able to locate it properly. Direct workaround will be to unzip .fasta.gz by yourself.
```
$ # At the directory where your data is located
$ gzip -d YOURFILE.fasta.gz # will generate YOURFILE.fasta
```
## Special Thanks & Disclaimer
1. This repository is a fork of https://github.com/SGSSGene/ImplementingSearch, which provided helpful skeleton code for this assignment. 
2. The dependencies stored in 'lib' subdirectory and all skeletal codes from above parent repository were accessed at the time of making this assignment (January 2024) and therefore will be possibly different from the updated version at the time of your view.
3. C++ Implementation was partly inspired by https://github.com/mbrner/ImplementingSearch. We strongly recommend reviewing this repository, as the author there employed many interesting approaches/coding style different from the code here.
4. The content of this repository was never endorsed by the authors of any of previosly mentioned repositories/libraries.
5. Our team name comes from https://www.youtube.com/watch?v=Oe3FG4EOgyU, which the author was endlessly listening to during the course of this assignment
