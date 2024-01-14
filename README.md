# ImplementingSearch
### Team "SchnappiDasKleineKrokodil"

Group 1: Antonio Alfaro de Prado, Eleanor Alspaugh, Hyunchang Oh

This is an exercise to demonstrate the power of Suffix-Array and FM-Index based searching.
Visit our Repository at https://github.com/HyunchangOh/ImplementingSearch

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

## Benchmark Comparison
For the query size of 1000,
Suffix array search took 343.96045565605164 seconds, among which 7.536531209945679 seconds were used for array construction.
Naive Search took 222.48520278930664 seconds.

## Discussion
It appears Python's native function 'find' is more time-efficient than our implementation of Search with suffix array.