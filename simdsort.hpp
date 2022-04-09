//
//
// simdsort.hpp
//
// 
// Copyright 2020- J. Makino all rights resesrved.
//
// sort with avx2 (etc)
#pragma once
#include <unistd.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#if defined(AVX2) || defined(AVX512) 
#include <x86intrin.h>
#endif

#if defined(AVX2) || defined(AVX512)  || defined(SVE)
#define SIMDSORTLIB_USE_SIMD
#endif

#if defined(SVE) 
#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>
#endif /* __ARM_FEATURE_SVE */
#endif

namespace SIMDSortLib{
#ifdef AVX2
    void make_permute_and_popcount_table_for_avx2(__m256i * stable,
						  int* ptable,
						  bool upper)
    {
	// lower 4 bits: 1 if smaller than pivot
	// upper 4 bits: 1 if larger than pivot
    
	for (int i=0;i<16;i++){
	    int count=0;
	    int mask=0;
	    int32_t * p = (int32_t*)(stable+i);
	    for(int k=0;k<4;k++){
		if ((1<<k)&i){
		    int index = count*2;;
		    if (upper){
			index = 6-index;
		    }
		    p[index] = k*2;
		    p[index+1] = k*2+1;
		    count++;
		}
	    }
	    ptable[i]=count;
	    //	fprintf(stderr, "i=%x count=%x ",i, count);
	    //	for (int ii=0;ii<8;ii++){
	    //	fprintf(stderr, " %x",p[ii]);
	    //	}
	    //	fprintf(stderr, "\n");
	}
    }
#endif
#ifdef AVX512
    void make_permute_and_popcount_table_for_avx512(__m512i * stable,
						    int* ptable,
						    bool upper)
    {
	// lower 4 bits: 1 if smaller than pivot
	// upper 4 bits: 1 if larger than pivot
    
	for (int i=0;i<256;i++){
	    int count=0;
	    int mask=0;
	    int64_t * p = (int64_t*)(stable+i);
	    for(int k=0;k<8;k++){
		if ((1<<k)&i){
		    int index = count;;
		    if (upper){
			index = 7-index;
		    }
		    p[index] = k;
		    count++;
		}
	    }
	    ptable[i]=count;
	    //fprintf(stderr, "i=%x count=%x ",i, count);
	    //	for (int ii=0;ii<8;ii++){
	    //	fprintf(stderr, " %lx",p[ii]);
	    //		}
	    //	fprintf(stderr, "\n");
	}
    }

    void make_mask_table_for_avx512(__mmask8 * table)
    {
	//  4 bits: 1 if larger than pivot
    
	for (int i=0;i<256;i++){
	    int count=0;
	    for(int k=0;k<8;k++){
		if ((1<<k)&i){
		    count++;
		}
	    }
	    int mask=0;
	    for(int k=0;k<count;k++){
		mask |= 1<<(7-k);
	    }
	    table[i]=mask;
	}
    }

    static __m512i permute_table_lower_avx512[256];
    static __m512i permute_table_upper_avx512[256];
    static int popcount_table_upper_avx512[256];
    static int popcount_table_lower_avx512[256];
    static __mmask8 mask_table_avx512[256];


#endif

#ifdef AVX2
    void make_single_permute_and_popcount_table_for_avx2(__m256i * stable,
							 int* ptable)
    {
	// lower 4 bits: 1 if smaller than pivot
	// upper 4 bits: 1 if larger than pivot
    
	for (int i=0;i<16;i++){
	    int count=0;
	    int countl=0;
	    int mask=0;
	    int32_t * p = (int32_t*)(stable+i);
	    for(int k=0;k<4;k++){
		if ((1<<k)&i){
		    int index = 6-count*2;
		    p[index] = k*2;
		    p[index+1] = k*2+1;
		    count++;
		}else{
		    int index = countl*2;
		    p[index] = k*2;
		    p[index+1] = k*2+1;
		    countl++;
		}
		
	    }
	    ptable[i]=count;
	    //	fprintf(stderr, "i=%x count=%x ",i, count);
	    //	for (int ii=0;ii<8;ii++){
	    //	    fprintf(stderr, " %x",p[ii]);
	    //	}
	    //	fprintf(stderr, "\n");
	}
    }


    void make_mask_table_for_avx2(__m256i * table)
    {
	//  4 bits: 1 if larger than pivot
    
	for (int i=0;i<16;i++){
	    int count=0;
	    for(int k=0;k<4;k++){
		if ((1<<k)&i){
		    count++;
		}
	    }
	    int64_t * p = (int64_t *) (table + i);
	    for(int k=0;k<4;k++){p[k]=0;}
	    for(int k=0;k<count;k++){p[3-k]=-1;}
	}
    }



    static __m256i permute_table_lower[16];
    static __m256i permute_table_upper[16];
    static int popcount_table_upper[16];
    static int popcount_table_lower[16];
    static __m256i mask_table[16];

    void dump_tables(__m256i * stable,
		     int* ptable)

    {
	for(int i=0;i<16;i++){
	    fprintf(stderr, "i=%x %x", i, ptable[i]);
	    int32_t * p = (int32_t*)(stable+i);
	    for (int ii=0;ii<8;ii++){
		fprintf(stderr, " %x",p[ii]);
	    }
	    fprintf(stderr, "\n");
	}
    }

    void dump_sort_table()
    {
	dump_tables(permute_table_upper, popcount_table_upper);
	dump_tables(permute_table_lower, popcount_table_lower);
    }

#endif
    void init_sort_table()
    {
#ifdef AVX2    
	make_permute_and_popcount_table_for_avx2(permute_table_upper,
						 popcount_table_upper,
						 true);
	make_permute_and_popcount_table_for_avx2(permute_table_lower,
						 popcount_table_lower,
						 false);
	make_mask_table_for_avx2(mask_table);
#endif    
#ifdef AVX512    
	make_permute_and_popcount_table_for_avx512(permute_table_upper_avx512,
						   popcount_table_upper_avx512,
						   true);
	make_permute_and_popcount_table_for_avx512(permute_table_lower_avx512,
						   popcount_table_lower_avx512,
						   false);
	make_mask_table_for_avx512(mask_table_avx512);
#endif    
    }
	


    inline double GetWtime()
    {
	struct timespec ts;
	if (clock_gettime(CLOCK_MONOTONIC,&ts )){
	    printf("GetWtime Failed\n");
	}
	return ((double) ts.tv_sec)+ ts.tv_nsec*1e-9;
    }

    static double time0;
    void init_timer()
    {
	time0=GetWtime();
    }
    void print_dt(std::string s)
    {
	double time=GetWtime();
	printf("%s %g\n",&(s[0]), time-time0);
	time0=GetWtime();
    }

    int compare_int64_t(const void *a, const void *b)
    {
	int64_t ua = *((int64_t*)a);
	int64_t ub = *((int64_t*)b);
	if (ua<ub){
	    return -1;
	}else if(ua==ub){
	    return 0;
	}else{
	    return 1;
	}
    }



    void sort_int64_array( int64_t * r, int lo, int up )
    {
	int i, j;
	int64_t tempr;
	//    fprintf(stderr, "called with %d %d\n", lo, up);
	while ( up>lo ) {
	    i = lo;
	    j = up;
	    tempr = r[lo];
	    /*** Split data in two ***/
	    while ( i<j ) {
		for ( ; r[j]> tempr; j-- );
		for ( r[i]=r[j]; i<j && r[i]<=tempr; i++ );
		r[j] = r[i];
	    }
	    r[i] = tempr;
	    /*** Sort recursively, the smallest first ***/
	    if ( i-lo < up-i ) {
		sort_int64_array(r,lo,i-1);  lo = i+1;
	    }else{
		sort_int64_array(r,i+1,up);  up = i-1;
	    }
	}
    }
    void sort_int64_array2( int64_t * r, int lo, int up )
    {
	if (up-lo<1) return;
	//    fprintf(stderr, "called with %d %d\n", lo, up);
	int i, j;
	int64_t tempr;
	i = lo;
	j = up;
	tempr = r[lo];
	/*** Split data in two ***/
	while ( i<j ) {
	    for ( ; r[j]> tempr; j-- );
	    for ( r[i]=r[j]; i<j && r[i]<=tempr; i++ );
	    r[j] = r[i];
	}
	r[i] = tempr;
	sort_int64_array2(r,lo,i-1); 
	sort_int64_array2(r,i+1,up); 
    }

    void sort_int64( int64_t * r, int n)
    {
	sort_int64_array(r, 0, n-1);
    }

#define SIMD_WIDTH 4

    int block_partition(int64_t* data, int64_t pibot, int n)
    {
	int64_t work[n];
#pragma omp simd
	//#pragma  GCC ivdep
	for(int i=0;i<n;i++){
	    work[i]=data[i];
	}
	int l= -1;
	int h=n;
	int i;
	for(i=0;i<n;i++){
	    if(work[i]< pibot){
		l++;
		data[l]=work[i];
	    }else if(work[i]> pibot){
		h--;
		data[h]=work[i];
	    }
	}
	for(i=l+1; i<h; i++){
	    data[i]=pibot;
	}

	return i-1;
    }

#ifdef AVX2
    typedef union m256di{
	__m256d d;
	__m256i i;
	__m256 f;
    }M256DI, *PM256DI;

#include "bitonic8.h"
    void dump256(__m256i* pdata,
		 char * s)
    {
	fprintf(stderr, "%s ", s);
	int64_t * p = (int64_t*) pdata;
	for(int i=0;i<4; i++){
	    fprintf(stderr, " %ld", p[i]);
	}
	fprintf(stderr, "\n");
    }


    int simd_partition_avx2(int64_t* data, int64_t pivot, int n, int*hi)
    {
	int n4 = (n+3)/4;
	int nb = n4*4;
	__m256i  work256[n4];
	int64_t *  work = (int64_t*) work256;
	__m256d* src = (__m256d*)data;
	__m256d* dest = (__m256d*)work;
	//    fprintf(stderr, "addresses = %lx %lx\n",(int64_t) src, (int64_t) dest);
	for(int i=0;i<n4;i++){
	    dest[i]=_mm256_loadu_pd((double*)(src+i));
	}
	int l= -1;
	int h=n;
	int i;
	int n4l = n/4;
	__m256i* pwork = (__m256i*)work;
	__m256i pivotv = _mm256_broadcastq_epi64(*((__m128i*)(&pivot)));
   
	for(int ii=0;ii<n4l;ii++){
#if 0
	    int i4 = ii<<2;
	    for(i=i4; i<i4+4; i++){
		if(work[i]> pivot){
		    h--;
		    data[h]=work[i];
		}else if(work[i]< pivot){
		    l++;
		    data[l]=work[i];
		}
	    }
#else
	    register union m256di u, u2; 
	    u.i =  _mm256_cmpgt_epi64(pivotv, pwork[ii]);
	    int maskl = _mm256_movemask_pd(u.d);
	    u.i =  _mm256_cmpgt_epi64( pwork[ii], pivotv);
	    int masku = _mm256_movemask_pd(u.d);
	    register union m256di lower, upper;
	    int dl = popcount_table_upper[maskl];
	    int dh = popcount_table_upper[masku];
	    //		int dl = _mm_popcnt_u32(maskl);
	    //		int dh = _mm_popcnt_u32(masku);
	    lower.f =_mm256_permutevar8x32_ps(*((__m256*)(pwork+ii)),
					      *((__m256i*)(permute_table_lower
							   +maskl)));
	    _mm256_storeu_pd((double*)(data+l+1), lower.d);
	    upper.f =_mm256_permutevar8x32_ps(*((__m256*)(pwork+ii)),
					      *((__m256i*)(permute_table_upper
							   +masku)));
	    _mm256_maskstore_pd((double*)(data+h-4), mask_table[masku],
				upper.d);
	    l+=dl;
	    h-=dh; 
#endif	
	}
	for(i=n4l*4;i<n;i++){
	    if(work[i]< pivot){
		l++;
		data[l]=work[i];
	    }else if(work[i]> pivot){
		h--;
		data[h]=work[i];
	    }
	}
	for(i=l+1; i<h; i++){
	    data[i]=pivot;
	}
	*hi=h;
	return l+1;
    }

    int simd_partition_2part_avx2(int64_t* data, int64_t pivot, int n)
    {
	int n4 = (n+3)/4;
	int nb = n4*4;
	__m256i  work256[n4];
	int64_t *  work = (int64_t*) work256;
	__m256d* src = (__m256d*)(data+1);
	__m256d* dest = (__m256d*)work;
	//    fprintf(stderr, "addresses = %lx %lx\n",(int64_t) src, (int64_t) dest);
	for(int i=0;i<n4;i++){
	    dest[i]=_mm256_loadu_pd((double*)(src+i));
	}
	int l= -1;
	int h=n;
	int i;
	int n4l = (n-1)/4;
	__m256i* pwork = (__m256i*)work;
	__m256i pivotv = _mm256_broadcastq_epi64(*((__m128i*)(&pivot)));
   
	for(int ii=0;ii<n4l;ii++){
#if 0
	    int i4 = ii<<2;
	    for(i=i4; i<i4+4; i++){
		if(work[i]> pivot){
		    h--;
		    data[h]=work[i];
		}else{
		    l++;
		    data[l]=work[i];
		}
	    }
#else
	    register union m256di u, u2; 
	    u.i =  _mm256_cmpgt_epi64( pwork[ii], pivotv);
	    int masku = _mm256_movemask_pd(u.d);
	    register union m256di lower, upper;
	    int dh = popcount_table_upper[masku];
	    int dl = 4-dh;
	    upper.f =_mm256_permutevar8x32_ps(*((__m256*)(pwork+ii)),
					      *((__m256i*)(permute_table_upper
							   +masku)));
	    _mm256_storeu_pd((double*)(data+l+1), upper.d);
	    _mm256_maskstore_pd((double*)(data+h-4), mask_table[masku],
				upper.d);
	    l+=dl;
	    h-=dh; 
#endif	
	}
	for(i=n4l*4;i<n-1;i++){
	    if(work[i]> pivot){
		h--;
		data[h]=work[i];
	    }else{
		l++;
		data[l]=work[i];
	    }
	}
	data[l+1]=pivot;
    
	return l+1;
    }

#endif


#ifdef AVX512
#ifndef AVX512_TYPE
    typedef union m512di{
	__m512d d;
	__m512i i;
	__m512 f;
    }M512DI,*PM512DI ;
#endif

#include "bitonic16.h"
#include "bitonic16kv.h"

    void dump512(__m512i* pdata,
		 char * s)
    {
	fprintf(stderr, "%s ", s);
	int64_t * p = (int64_t*) pdata;
	for(int i=0;i<8; i++){
	    fprintf(stderr, " %lx", p[i]);
	}
	fprintf(stderr, "\n");
    }

    int simd_partition_avx512(int64_t* data, int64_t pivot, int n, int*hi)
    {
	//    if (n < 8){
	//	return simd_partition_avx2(data,  pivot,  n);
	//    }
	int n8 = (n+7)/8;
	int nb = n8*8;
	__m512i  work256[n8];
	int64_t *  work = (int64_t*) work256;
	__m512d* src = (__m512d*)data;
	__m512d* dest = (__m512d*)work;
	for(int i=0;i<n8;i++){
	    dest[i]=_mm512_loadu_pd((double*)(src+i));
	}
	int l= -1;
	int h=n;
	int i;
	int n8l = n/8;
	__m512i* pwork = (__m512i*)work;
	__m512i pivotv = _mm512_broadcastq_epi64(*((__m128i*)(&pivot)));
   
	for(int ii=0;ii<n8l;ii++){
#if 0
	    int i8 = ii<<3;
	    for(i=i8; i<i8+8; i++){
		if(work[i]> pivot){
		    h--;
		    dat<a[h]=work[i];
		}else if(work[i]< pivot){
		    l++;
		    data[l]=work[i];
		}
	    }
#else
	    __mmask8 maskl =  _mm512_cmpgt_epi64_mask(pivotv, pwork[ii]);
	    __mmask8 masku =  _mm512_cmpgt_epi64_mask( pwork[ii], pivotv);
	    int dl = _mm_popcnt_u32(maskl);
	    int dh = _mm_popcnt_u32(masku);
	    _mm512_mask_compressstoreu_pd((double*)(data+l+1), maskl, 
					  *((__m512d*)(pwork+ii)));
	    _mm512_mask_compressstoreu_pd((double*)(data+h-dh), masku, 
					  *((__m512d*)(pwork+ii)));
	    l+=dl;
	    h-=dh; 
#endif	
	}
	int nremain = n -n8l*8;
	if (nremain > 0){
	    int ii = n8l;
	    __mmask8 maskr =(__mmask8)  ((1<<nremain)-1);
	    __mmask8 maskl =  _mm512_mask_cmpgt_epi64_mask(maskr,pivotv, pwork[ii]);
	    __mmask8 masku =  _mm512_mask_cmpgt_epi64_mask(maskr, pwork[ii], pivotv);
	    int dl = _mm_popcnt_u32(maskl);
	    int dh = _mm_popcnt_u32(masku);
	    _mm512_mask_compressstoreu_pd((double*)(data+l+1), maskl, 
					  *((__m512d*)(pwork+ii)));
	    _mm512_mask_compressstoreu_pd((double*)(data+h-dh), masku, 
					  *((__m512d*)(pwork+ii)));
	    l+=dl;
	    h-=dh; 
	}
#if 0    
	for(i=n8l*8;i<n;i++){
	    if(work[i]< pivot){
		l++;
		data[l]=work[i];
	    }else if(work[i]> pivot){
		h--;
		data[h]=work[i];
	    }
	}
#endif    
	for(i=l+1; i<h; i++){
	    data[i]=pivot;
	}
	*hi = h;
	return l+1;
    }

    int simd_partition_avx512(uint64_t* vhi,
			      uint64_t* vlo,
			      uint64_t* index,
			      uint64_t pivothi,
			      uint64_t pivotlo,
			      uint64_t pivotindex,
			      int n,
			      int*hi)
    {
	int n8 = (n+7)/8;
	int nb = n8*8;
	__m512i  work512hi[n8];
	__m512i  work512lo[n8];
	__m512i  work512index[n8];
	uint64_t *  workhi = (uint64_t*) work512hi;
	uint64_t *  worklo = (uint64_t*) work512lo;
	uint64_t *  workindex = (uint64_t*) work512index;
	__m512d* srchi = (__m512d*)vhi;
	__m512d* desthi = (__m512d*)workhi;
	__m512d* srclo = (__m512d*)vlo;
	__m512d* destlo = (__m512d*)worklo;
	__m512d* srcindex = (__m512d*)index;
	__m512d* destindex = (__m512d*)workindex;
	for(int i=0;i<n8;i++){
	    desthi[i]=_mm512_loadu_pd((double*)(srchi+i));
	    destlo[i]=_mm512_loadu_pd((double*)(srclo+i));
	    destindex[i]=_mm512_loadu_pd((double*)(srcindex+i));
	}
	int l= -1;
	int h=n;
	int i;
	int n8l = n/8;
	__m512i* pworkhi = (__m512i*)workhi;
	__m512i* pworklo = (__m512i*)worklo;
	__m512i* pworkindex = (__m512i*)workindex;
	__m512i pivotvhi = _mm512_broadcastq_epi64(*((__m128i*)(&pivothi)));
	__m512i pivotvlo = _mm512_broadcastq_epi64(*((__m128i*)(&pivotlo)));
	__m512i pivotvindex = _mm512_broadcastq_epi64(*((__m128i*)(&pivotindex)));
   
	for(int ii=0;ii<n8l;ii++){
   
	    __mmask8 maskhieq =  _mm512_cmpeq_epu64_mask(pivotvhi, pworkhi[ii]);
	    __mmask8 maskl =  _mm512_cmpgt_epu64_mask(pivotvhi, pworkhi[ii])|
		(maskhieq &  _mm512_cmpgt_epu64_mask(pivotvlo, pworklo[ii]));
	    __mmask8 masku =  _mm512_cmplt_epu64_mask(pivotvhi, pworkhi[ii])|
		(maskhieq &  _mm512_cmplt_epu64_mask(pivotvlo, pworklo[ii]));
	    
	    int dl = _mm_popcnt_u32(maskl);
	    int dh = _mm_popcnt_u32(masku);
	    _mm512_mask_compressstoreu_pd((double*)(vhi+l+1), maskl, 
					  *((__m512d*)(pworkhi+ii)));
	    _mm512_mask_compressstoreu_pd((double*)(vlo+l+1), maskl, 
					  *((__m512d*)(pworklo+ii)));
	    _mm512_mask_compressstoreu_pd((double*)(index+l+1), maskl, 
					  *((__m512d*)(pworkindex+ii)));
	    _mm512_mask_compressstoreu_pd((double*)(vhi+h-dh), masku, 
					  *((__m512d*)(pworkhi+ii)));
	    _mm512_mask_compressstoreu_pd((double*)(vlo+h-dh), masku, 
					  *((__m512d*)(pworklo+ii)));
	    _mm512_mask_compressstoreu_pd((double*)(index+h-dh), masku, 
					  *((__m512d*)(pworkindex+ii)));
	    l+=dl;
	    h-=dh; 
	}
	int nremain = n -n8l*8;
	if (nremain > 0){
	    int ii = n8l;
	    __mmask8 maskr =(__mmask8)  ((1<<nremain)-1);
	    __mmask8 maskhieq =  _mm512_cmpeq_epu64_mask(pivotvhi, pworkhi[ii]);
	    __mmask8 maskl =  _mm512_cmpgt_epu64_mask(pivotvhi, pworkhi[ii])|
		(maskhieq &  _mm512_cmpgt_epu64_mask(pivotvlo, pworklo[ii]));
	    __mmask8 masku =  _mm512_cmplt_epu64_mask(pivotvhi, pworkhi[ii])|
		(maskhieq &  _mm512_cmplt_epu64_mask(pivotvlo, pworklo[ii]));
	    maskl &= maskr;
	    masku &= maskr;
	    int dl = _mm_popcnt_u32(maskl);
	    int dh = _mm_popcnt_u32(masku);
	    _mm512_mask_compressstoreu_pd((double*)(vhi+l+1), maskl, 
					  *((__m512d*)(pworkhi+ii)));
	    _mm512_mask_compressstoreu_pd((double*)(vlo+l+1), maskl, 
					  *((__m512d*)(pworklo+ii)));
	    _mm512_mask_compressstoreu_pd((double*)(index+l+1), maskl, 
					  *((__m512d*)(pworkindex+ii)));
	    _mm512_mask_compressstoreu_pd((double*)(vhi+h-dh), masku, 
					  *((__m512d*)(pworkhi+ii)));
	    _mm512_mask_compressstoreu_pd((double*)(vlo+h-dh), masku, 
					  *((__m512d*)(pworklo+ii)));
	    _mm512_mask_compressstoreu_pd((double*)(index+h-dh), masku, 
					  *((__m512d*)(pworkindex+ii)));
	    l+=dl;
	    h-=dh; 
	}
	for(i=l+1; i<h; i++){
	    vhi[i]=pivothi;
	    vlo[i]=pivotlo;
	    index[i]=pivotindex;
	}
	*hi = h;
	return l+1;
    }


#endif

#ifdef SVE
#include "bitonic16sve.h"
#include "bitonic16svekv.h"

    void dumpsve(svint64_t v,
		 char *s)
    {
	fprintf(stderr, "%s ", s);
	int64_t x[8];
	svst1_s64(svptrue_b64(),x, v);
	for(int i=0;i<8;i++){
	    fprintf(stderr, " %lx", x[i]);
	}
	fprintf(stderr, "\n");
    
    }

    void dump512(int64_t * pdata,
		 char * s)
    {
	fprintf(stderr, "%s ", s);
	int64_t * p = (int64_t*) pdata;
	for(int i=0;i<8; i++){
	    fprintf(stderr, " %lx", p[i]);
	}
	fprintf(stderr, "\n");
    }


    int simd_partition_sve(int64_t* data, int64_t pivot, int n, int* hi)
    {
	int nsve = svcntd();
	int n8 = (n+nsve)/nsve;
	int nb = n8*nsve;
	int64_t __attribute__ ((aligned(64)))  work[nb];
	svbool_t ptrue =svptrue_b64();
	for(int i=0;i<nb;i+=nsve){
	    svst1_s64(ptrue, work+i, svld1_s64(ptrue, data+i));
	}
	int l= -1;
	int h=n;
	int i;
	int n8l = n/nsve;
	svint64_t  pivotv = svdup_s64(pivot);
   
#if 0
	for(int i8=0;i8<n8l*nsve;i8+=nsve){
	    for(i=i8; i<i8+8; i++){
		if(work[i]> pivot){
		    h--;
		    dat<a[h]=work[i];
		}else if(work[i]< pivot){
		    l++;
		    data[l]=work[i];
		}
	    }
	}
	for(i=n8l*8;i<n;i++){
	    if(work[i]< pivot){
		l++;
		data[l]=work[i];
	    }else if(work[i]> pivot){
		h--;
		data[h]=work[i];
	    }
	}
#else
	for(int i8=0;i8<n;i8+=nsve){
	    svbool_t pmask= svwhilelt_b64(i8,n);
	    svint64_t val =svld1_s64(pmask, work+i8);
	    svbool_t maskl =  svcmpgt_s64(pmask, pivotv,val);
	    svbool_t masku =  svcmpgt_s64(pmask, val, pivotv);
	    int dl = svcntp_b64(pmask, maskl);
	    int dh = svcntp_b64(pmask, masku);
	    svst1_s64(pmask, data+l+1,  svcompact_s64(maskl, val));
	    svbool_t maskustore = svwhilelt_b64(0,dh);
	    svst1_s64(maskustore, data+h-dh,  svcompact_s64(masku, val));
	    l+=dl;
	    h-=dh; 
	}
#endif	
	for(i=l+1; i<h; i++){
	    data[i]=pivot;
	}
	
    
	*hi = h;
	return l+1;
    }

    int simd_partition_sve(uint64_t* vhi,
			      uint64_t* vlo,
			      uint64_t* index,
			      uint64_t pivothi,
			      uint64_t pivotlo,
			      uint64_t pivotindex,
			      int n,
			      int*hi)
    {
	int nsve = svcntd();
	int n8 = (n+nsve)/nsve;
	int nb = n8*nsve;
	uint64_t __attribute__ ((aligned(64)))  workhi[nb];
	uint64_t __attribute__ ((aligned(64)))  worklo[nb];
	uint64_t __attribute__ ((aligned(64)))  workindex[nb];
	svbool_t ptrue =svptrue_b64();
	for(int i=0;i<nb;i+=nsve){
	    svst1(ptrue, workhi+i, svld1(ptrue, vhi+i));
	    svst1(ptrue, worklo+i, svld1(ptrue, vlo+i));
	    svst1(ptrue, workindex+i, svld1(ptrue, index+i));
	}
	int l= -1;
	int h=n;
	int i;
	int n8l = n/nsve;
	svuint64_t  pivotvhi = svdup_u64(pivothi);
	svuint64_t  pivotvlo = svdup_u64(pivotlo);
	svuint64_t  pivotvindex = svdup_u64(pivotindex);
   
	for(int i8=0;i8<n;i8+=nsve){
	    svbool_t pmask= svwhilelt_b64(i8,n);
	    svuint64_t valhi =svld1(pmask, workhi+i8);
	    svuint64_t vallo =svld1(pmask, worklo+i8);
	    svuint64_t valindex =svld1(pmask, workindex+i8);
	    svbool_t masklhi =  svcmpgt(pmask, pivotvhi,valhi);
	    svbool_t maskuhi =  svcmpgt(pmask, valhi, pivotvhi);
	    svbool_t maskeqhi =  svcmpeq(pmask, valhi, pivotvhi);
	    svbool_t maskllo =  svcmpgt(pmask, pivotvlo,vallo);
	    svbool_t maskulo =  svcmpgt(pmask, vallo, pivotvlo);
	    svbool_t maskeqlo =  svcmpeq(pmask, vallo, pivotvlo);
	    svbool_t masklindex =  svcmpgt(pmask, pivotvindex,valindex);
	    svbool_t maskuindex =  svcmpgt(pmask, valindex, pivotvindex);
	    svbool_t maskl = svorr_z(pmask,masklhi,
				   svorr_z(pmask,
					 svand_z(pmask,maskeqhi, maskllo),
					 svand_z(pmask,maskeqhi,
					       svand_z(pmask, maskeqlo,
						     masklindex))));
	    svbool_t masku = svorr_z(pmask,maskuhi,
				   svorr_z(pmask,
					 svand_z(pmask,maskeqhi, maskulo),
					 svand_z(pmask,maskeqhi,
					       svand_z(pmask, maskeqlo,
						     maskuindex))));
	    
								     

	    int dl = svcntp_b64(pmask, maskl);
	    int dh = svcntp_b64(pmask, masku);
	    svst1(pmask, vhi+l+1,  svcompact(maskl, valhi));
	    svst1(pmask, vlo+l+1,  svcompact(maskl, vallo));
	    svst1(pmask, index+l+1,  svcompact(maskl, valindex));
	    svbool_t maskustore = svwhilelt_b64(0,dh);
	    svst1(maskustore, vhi+h-dh,  svcompact(masku, valhi));
	    svst1(maskustore, vlo+h-dh,  svcompact(masku, vallo));
	    svst1(maskustore, index+h-dh,  svcompact(masku, valindex));
	    l+=dl;
	    h-=dh; 
	}
	for(i=l+1; i<h; i++){
	    vhi[i]=pivothi;
	    vlo[i]=pivotlo;
	    index[i]=pivotindex;
	}
	
    
	*hi = h;
	return l+1;
    }

#endif

    int simd_partition(uint64_t* vhi,
		       uint64_t* vlo,
		       uint64_t* index,
		       int n, int* hi)
    {
	//	int nsve = svcntd();
	int nsve = 8;
	int n8 = (n+nsve)/nsve;
	int nb = n8*nsve;
	int64_t __attribute__ ((aligned(64)))  workhi[nb];
	int64_t __attribute__ ((aligned(64)))  worklo[nb];
	int64_t __attribute__ ((aligned(64)))  workindex[nb];
	for(int i=0;i<n;i++){
	    workhi[i]=vhi[i];
	    worklo[i]=vlo[i];
	    workindex[i]=index[i];
	}
	int l= -1;
	int h=n;
	int i;
	int n8l = n/nsve;
	uint64_t pivothi = workhi[0];
	uint64_t pivotlo = worklo[0];
	uint64_t pivotindex = workindex[0];
	for(i=0;i<n;i++){
	    if(workhi[i]< pivothi
	       || (workhi[i]== pivothi && worklo[i]< pivotlo )
	       || (workhi[i]== pivothi &&
		   worklo[i]== pivotlo && workindex[i]<pivotindex)){
		l++;
		vhi[l]=workhi[i];
		vlo[l]=worklo[i];
		index[l]=workindex[i];
	    }else if(workhi[i]> pivothi
		     || (workhi[i]== pivothi && worklo[i]> pivotlo )
		     || (workhi[i]== pivothi &&
			 worklo[i]== pivotlo && workindex[i]>pivotindex)){
		h--;
		vhi[h]=workhi[i];
		vlo[h]=worklo[i];
		index[h]=workindex[i];
	    }
	}
	for(i=l+1; i<h; i++){
	    vhi[i]=pivothi;
	    vlo[i]=pivotlo;
	    index[i]=pivotindex;
	}
	
    
	*hi = h;
	return l+1;
    }

    

    void dump_data(int64_t* data,
		   int n,
		   char* message)
    {
	printf("%s ", message);
	for(int i=0;i<n; i++){
	    printf(" %lx", data[i]);
	}
	printf("\n");
    }

    void simd_sort_recursivex( int64_t * r, int lo, int up ) 
    {
	//    fprintf(stderr, "called with %d %d\n", lo, up);
	if (up-lo<1) return;
#ifdef AVX2
	if (up-lo+1<=8){
	    //	fprintf(stderr, "call bitonic8 with  %d\n", up-lo+1);
	
	    bitonic8(r+lo, up-lo+1);
	    return;
	}
	    
#endif	    
#if defined(AVX512) || defined(SVE)
	if (up-lo+1<=16){
	    //	fprintf(stderr, "call bitonic16 with  %d\n", up-lo+1);
	
	    bitonic16(r+lo, up-lo+1);
	    return;
	}
	    
#endif	    
	int i, j, hi;
	int64_t tempr;
	//    dump_data( r+lo, up-lo+1, "before");
#ifdef AVX2    
	i=simd_partition_avx2(r+lo, r[lo], up-lo+1, &hi);
#endif
#ifdef AVX512
	i=simd_partition_avx512(r+lo, r[lo], up-lo+1, &hi);
#endif
#ifdef SVE
	i=simd_partition_sve(r+lo, r[lo], up-lo+1, &hi);
#endif
    
	//   dump_data( r+lo, up-lo+1, "after ");
	//    printf("i=%d\n", i);
	simd_sort_recursivex(r,lo,lo+i-1);  
	simd_sort_recursivex(r,lo+hi,up);  
    
    }


    void simd_sort_recursive(uint64_t* vhi,
			     uint64_t* vlo,
			     uint64_t* index,
			     int lo,
			     int up )
    {
	//    fprintf(stderr, "called with %d %d\n", lo, up);
	if (up-lo<1) return;
#ifdef AVX2XX
	if (up-lo+1<=8){
	    //	fprintf(stderr, "call bitonic8 with  %d\n", up-lo+1);
	
	    bitonic8(r+lo, up-lo+1);
	    return;
	}
	    
#endif	    
#if defined(AVX512) || defined(SVE)
	if (up-lo+1<=16){
	    //	fprintf(stderr, "call bitonic16 with  %d\n", up-lo+1);
	
	    bitonic16(vhi+lo,vlo+lo, index+lo,  up-lo+1);
	    return;
	}
	    
#endif	    
	int i, j, hi;
	int64_t tempr;
	//    dump_data( r+lo, up-lo+1, "before");
#ifdef AVX2XX    
	i=simd_partition_avx2(r+lo, r[lo], up-lo+1, &hi);
#endif
#ifdef AVX512
	i=simd_partition_avx512(vhi+lo, vlo+lo, index+lo,
				vhi[lo],vlo[lo], index[lo],
				up-lo+1, &hi);
#endif
#ifdef SVE
	i=simd_partition_sve(vhi+lo, vlo+lo, index+lo,
				vhi[lo],vlo[lo], index[lo],
				up-lo+1, &hi);

#endif
#ifndef SIMDSORTLIB_USE_SIMD	
	i=simd_partition(vhi+lo, vlo+lo, index+lo,   up-lo+1, &hi);
#endif    
	//   dump_data( r+lo, up-lo+1, "after ");
	//    printf("i=%d\n", i);
	simd_sort_recursive(vhi, vlo, index, lo,lo+i-1);  
	simd_sort_recursive(vhi, vlo, index, lo+hi,up);  
    
    }



 
    void simd_sort( int64_t * r, int n)
    {
	static int initialized=0;
#ifdef SIMDSORTLIB_USE_SIMD	
	if (initialized==0){
	    init_sort_table();
#ifdef AVX2
	    initialize_lsmask();
#endif	    
	    initialized=1;
	}
#endif	    
	simd_sort_recursivex(r, 0, n-1);
    }

    void simd_sort( uint64_t* hi,uint64_t* lo, uint64_t* index,  int n)
    {
	static int initialized=0;
	if (initialized==0){
	    init_sort_table();
	    initialized=1;
	}
	simd_sort_recursive(hi, lo, index,  0, n-1);
    }
    
	
}

 
