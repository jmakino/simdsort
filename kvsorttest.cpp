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
#include<string>
#include<map>
#include<random>
#include<omp.h>
#include<stdlib.h>


#include "simdsort.hpp"

class KeyValuePair{
public:
    uint64_t key;
    uint64_t hi;
    uint64_t lo;
    inline uint64_t get_hi_key(){return hi;}
    inline uint64_t get_lo_key(){return lo;}
    inline uint64_t get_index(){return key;}
    inline void  set_hi_key(uint64_t val){hi=val;}
    inline void  set_lo_key(uint64_t val){lo=val;}
    inline void  set_index(uint64_t val){key=val;}
};


template <class T>
auto has_get_hi_key(T& a) -> decltype(a.get_hi_key(), bool())
{
  return true;
}

auto has_get_hi_key(...) -> bool
{
  return false;
}


template<class T>
bool operator<(T& a, T& b)
{
    if (a.get_hi_key() < b.get_hi_key() ||  ((a.get_hi_key() == b.get_hi_key()) && (a.get_lo_key() < b.get_lo_key()))){
	return true;
    }else{
	return false;
    }
}

template<class T>
bool operator==(T& a, T& b)
{
    return  (a.get_hi_key() == b.get_hi_key() &&  a.get_lo_key() ==  b.get_lo_key() &&  a.key ==  b.key);
}

template<class T>
bool operator!=(T& a, T& b)
{
    return  ! (a==b);
}

void dump_data(KeyValuePair * data,
	       int n,
	       std::string s)
{
    fprintf(stderr, "%s:\n", &(s[0]));
    for(auto i=0;i<n;i++){
	fprintf(stderr, "data[%d] = %lu %lu %lu\n",
		i, data[i].get_hi_key(), data[i].get_lo_key(),  data[i].key);
    }
}

    
template<class T>
bool operator>(T& t1, T& t2) { return t2 < t1; }
template<class T>
bool operator<=(T& t1, T& t2) { return !(t1 > t2); }
template<class T>
bool operator>=(T& t1, T& t2) { return !(t1 < t2); }
    
int main(int argc, char** argv)
{
    auto nstart= atoi(argv[1]);
    auto nend =nstart+1;
    if (argc > 2) nend= atoi(argv[2]);
    bool showtime = false;
    bool showdata = false;
    if (argc > 3) {
	showtime = (atoi(argv[3]) == 1);
	showdata = (atoi(argv[3]) == 2);
    }else if (argc <= 2){
	showtime = true;
    }
    
    KeyValuePair data[nend];
    KeyValuePair data0[nend];
    KeyValuePair data1[nend];

    
    srandom(1);
    for (int n=nstart; n<nend; n++){
	for (int i=0; i<n; i++){
	    data[i].hi=(random()<<32 | random()) &3; 
	    data[i].lo=(random()<<32 | random()) ; 
	    data[i].key=i;
	    data0[i]=data[i];
	    data1[i]=data[i];
	}
	//	for (int i=0; i<n; i++){
	//	    printf("i=%d   %ld\n",i, data[i]);
	//	}
	if (showdata) dump_data(data, n, "before sort");

	bool ok=true;
	SIMDSortLib::init_timer();
	std::sort(data, data+n); 
	if (showtime)SIMDSortLib::print_dt("C++ std::sort:  ");
	SIMDSortLib::init_timer();
	SIMDSortLib::sort_kv_using_simdsort(data1, n); 
	if (showtime)SIMDSortLib::print_dt("blocksort: ");
	//	for (int i=0; i<n; i++){
	//	    printf("i=%d   %ld %ld\n",i, data[i], data1[i]);
	//	}

	if (showdata) dump_data(data, n, "after std sort");
	if (showdata) dump_data(data1, n, "after sort");
	for (int i=0; i<n; i++){
	    if (data[i] != data1[i]){
		ok=false;
		printf("ERROR at:%d  %lu:%lu   %lu:%lu \n",i,
		       data[i].hi, data[i].lo,
		       data1[i].hi, data1[i].lo);
		exit(-1);
	    }
	}
	if (ok){
	    printf("sort passed for n=%d\n", n);
	}else{
	    printf("Parallel sort failed!\n");
	}

    }
    printf("data has hi key=%d\n", has_get_hi_key(data[0]));
    int n;
    printf("int has hi key=%d\n", has_get_hi_key(n));
    exit(0);
    
}

