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

## Benchmark Comparison
For the query size of 1000,
Suffix array search took 343.96045565605164 seconds, among which 7.536531209945679 seconds were used for array construction.
Naive Search took 222.48520278930664 seconds.

## Time and Memory Usage Comparison
For query size of 10,100,1000, the memory usage and computation time of each search method was compared, using python's native time library for computation time and memory-profiler for memory usage.
For detailed use of memory-profiler library, visit https://pypi.org/project/memory-profiler/


### Memory Usage
![image](https://github.com/HyunchangOh/ImplementingSearch/assets/42934606/55784363-5bbe-4971-ab44-d01480e5e11d)
Naive search does not generate any additional data structures, so its memory use is the lowest, while Suffix Array uses significantly more memory than FM index.

### Computation Time
![image](https://github.com/HyunchangOh/ImplementingSearch/assets/42934606/a6964386-2e58-411c-a589-f50922a9c890)
Of the three algorithms tested between assignments 1 and 2, the FMindex search appears to be the most efficient. One reason for this is its average-case linear time complexity. The FMindex search algorithm also uses a log-based time complexity for encoding of the test while suffix array is not log-based; this makes it faster than suffix array for individual searches. It is based on the Burrows Wheeler Transform (BWT) and suffix array. The FMindex search also has lower time requirements than naive search. The naive search does not require the same pre-processing that FMindex and suffix array require but this pre-processing step reduces downstream processing time during the search.

# C++ Program
## Dependencies
1. This program runs with Seqan3 and libsufsort, which is already installed in 'lib' subdirectory.
2. The dependencies, especially seqan3, require specific versions of compilers (and possibly others). Check troubleshooting if you run into errors upon running 'cmake'.

## Implementation Description
### Naive Search

### Suffix Search

### Fm Index Search

### Pigeonhole Search

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

## Troubleshooting
### Pushing a Large File to a Github Repo
If pushing to a github repository returns:
fatal: the remote end hung up unexpectedly

Try (but it is not recommended to push a large file like .fasta file directly to a github repo):
git config http.postBuffer 524288000

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