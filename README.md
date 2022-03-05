# simdsort

A fast quicksort function with AVX2 (currently) support

you can try this with

```
  make testsort
```
The last part of the  output would look like
```
./simdsorttest 8001
qsort:  0.000525483
sort:   0.000371819
blocksort:  0.000224667
sort passed for n=8001
```
The numbers are elapsed time to sort 8001 64-bit integer numbers,
with system qsort(3), non-simd quiksort and simd quicksort.
This result is ontained  on Intel® Core™ i7-1065G7 with gcc version
7.5.0 (Ubuntu 7.5.0-3ubuntu1~18.04).  



## Usage as header-only library

```
   #include "PATH/OF/THIS/FILE/simdsort.h"
```



Provide:
```
void simd_sort_int64( int64_t * r, int n);
```


r is the pointer to the array of int64_t to be sorted
n is the number of elements

## Note

This source is in pure C and no care for name collision is done.

## Algorithm

Essentially the same as in
Fast Quicksort Implementation Using AVX Instructions,
Gueron and  Krasnov 2015, 
10.1093/comjnl/bxv063
