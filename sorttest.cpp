#include<cmath>
#include<vector>
#include<functional>
#include<algorithm>
#include<exception>
#include<stdexcept>
#include<cassert>
#include<typeinfo>
#include<cstdio>
#include<cstring>
#include<map>
#include<random>
#include<omp.h>
#include<stdlib.h>


#include "simdsort.hpp"

int main(int argc, char** argv)
{
    auto nstart= atoi(argv[1]);
    auto nend =nstart+1;
    //    using namespace SampleSortLib;
    // 一様実数分布
    // [-1.0, 1.0)の値の範囲で、等確率に実数を生成する
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
	SIMDSortLib::init_timer();
	std::sort(data, data+n); 
	if (showtime)SIMDSortLib::print_dt("C++ std::sort:  ");
	SIMDSortLib::init_timer();
	SIMDSortLib::simd_sort(data1, n); 
	if (showtime)SIMDSortLib::print_dt("blocksort: ");
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

