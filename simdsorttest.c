#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
//#include <x86intrin.h>

#include "simdsort.h"

void sort_int64_c(int64_t * a, int length);//defined in sort.cpp


int main(int argc, char** argv)
{
    //    init_sort_table();
    //    dump_sort_table();
    int nstart= atoi(argv[1]);
    int nend =nstart+1;
    if (argc > 2) nend= atoi(argv[2]);
    bool showtime = false;
    if (argc > 3) {
	showtime = (atoi(argv[3]) == 1);
    }else if (argc <= 2){
	showtime = true;
    }
    
    int64_t data[nend];
    int64_t data0[nend];
    int64_t data1[nend];

    
    srandom(1);
    for (int n=nstart; n<nend; n++){
	for (int i=0; i<n; i++){
	    data[i]=(random()<<32 | random()) ; 
	    data0[i]=data[i];
	    data1[i]=data[i];
	}
	//	for (int i=0; i<n; i++){
	//	    printf("i=%d   %ld\n",i, data[i]);
	//	}

	bool ok=true;
	init_timer();
	qsort(data, n, sizeof(int64_t), compare_int64_t);
	if (showtime)print_dt("qsort: ");
	init_timer();
	sort_int64_c(data0, n); 
	if (showtime)print_dt("C++ std::sort:  ");
	init_timer();
	simd_sort_int64(data1, n); 
	if (showtime)print_dt("blocksort: ");
	//	for (int i=0; i<n; i++){
	//	    printf("i=%d   %ld %ld\n",i, data[i], data1[i]);
	//	}

	for (int i=0; i<n; i++){
	    if (data[i] != data1[i]){
		ok=false;
		printf("ERROR at:%d  %ld %ld\n",i, data[i], data1[i]);
		exit(-1);
	    }
	}
	if (ok){
	    printf("sort passed for n=%d\n", n);
	}else{
	    printf("Parallel sort failed!\n");
	}
    }
    exit(0);
    
}

