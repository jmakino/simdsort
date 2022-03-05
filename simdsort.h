//
//
// simdsort.h
//
// 
// Copyright 2020- J. Makino all rights resesrved.
//
// sort with avx2 (etc)

#include <stdbool.h>
#include <assert.h>
#include <x86intrin.h>
#include <time.h>

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

void init_sort_table()
{
    make_permute_and_popcount_table_for_avx2(permute_table_upper,
					     popcount_table_upper,
					     true);
    make_permute_and_popcount_table_for_avx2(permute_table_lower,
					     popcount_table_lower,
					     false);
    make_mask_table_for_avx2(mask_table);
}
	
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
void print_dt(char* s)
{
    double time=GetWtime();
    printf("%s %g\n", s, time-time0);
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

union m256di{
    __m256d d;
    __m256i i;
    __m256 f;
};

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
    
	    

int simd_partition_avx2(int64_t* data, int64_t pivot, int n)
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
	    if(work[i]< pivot){
		l++;
		data[l]=work[i];
	    }else if(work[i]> pivot){
		h--;
		data[h]=work[i];
	    }
	} 
#else
	union m256di u, u2; 
	u.i =  _mm256_cmpgt_epi64(pivotv, pwork[ii]);
	int maskl = _mm256_movemask_pd(u.d);
	u.i =  _mm256_cmpgt_epi64( pwork[ii], pivotv);
	int masku = _mm256_movemask_pd(u.d);
	union m256di lower, upper;
	int dl = popcount_table_upper[maskl];
	int dh = popcount_table_upper[masku];
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
    
    return i-1;
}

void dump_data(int64_t* data,
	       int n,
	       char* message)
{
    printf("%s ", message);
    for(int i=0;i<n; i++){
	printf(" %ld", data[i]);
    }
    printf("\n");
}


void simd_sort_int64_array( int64_t * r, int lo, int up )
{
    //    printf("called with %d %d\n", lo, up);
    if (up-lo<1) return;
    int i, j;
    int64_t tempr;
    //    dump_data( r+lo, up-lo+1, "before");
    i=simd_partition_avx2(r+lo, r[lo], up-lo+1);
    //    dump_data( r+lo, up-lo+1, "after ");
    //    printf("i=%d\n", i);
    simd_sort_int64_array(r,lo,lo+i-1);  
    simd_sort_int64_array(r,lo+i+1,up);  
}

void simd_sort_int64( int64_t * r, int n)
{
    static int initialized=0;
    if (initialized==0){
	init_sort_table();
	initialized=1;
    }
    simd_sort_int64_array(r, 0, n-1);
}

    
	
