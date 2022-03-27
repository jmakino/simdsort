//
// bitonic8.h
#pragma once
typedef uint64_t INT64T;

void permute(INT64T *dest,
	     INT64T *src,
	     INT64T *dest_table,
	     int n)
{
    for(int i=0;i<n;i++){
	dest[dest_table[i]]=src[i];
    }
}
void copy(INT64T *dest,
	  INT64T *src,
	      int n)
{
    for(int i=0;i<n;i++){
	dest[i]=src[i];
    }
}
void dump(INT64T * data, char * s)
{
    fprintf(stderr, "%s", s);
    for(int i=0;i<8;i++){
	fprintf(stderr, " %lu", data[i]);
    }
    fprintf(stderr, "\n");
}

void dumpavx2(M256DI a,
	      M256DI b,
	      char *s)
{
    INT64T data[8];
    _mm256_store_pd((double*)data, a.d);
    _mm256_store_pd((double*)(data+4), b.d);
    dump(data, s);
}

void compare(INT64T * data,
	     int n)
{
    int nhalf=n/2;
    for(int i=0;i<nhalf; i++){
	INT64T a, b;
	a=data[i];
	b=data[i+nhalf];
	data[i] = (a<b)? a:b;
	data[i+nhalf] = (a>=b)? a:b;
    }
}
void bitonic8_basic(INT64T * data)
{
    INT64T work[8];

    copy(work,data,8);
    
    dump(work, "after first permute");
    compare(work,8);
    dump(work, "after compare");
    INT64T tab2[]={0,1,2,3,5,4,7,6};
    permute(data,work, tab2,8);  //1
    dump(data, "after second permute");
    compare(data,8);
    dump(data, "after 2nd compare");
    INT64T tab2a[]={0,4,2,6,5,1,7,3};
    permute(work,data, tab2a,8);  //2
    dump(work, "after 3rd permute");
    compare(work,8);
    dump(work, "after 3rd comp");
    
    INT64T tab3[]={0,2,7,5,1,3,6,4};
    permute(data,work, tab3,8);  //3
    dump(data, "after 4th permute");
    compare(data,8);
    INT64T tab4[]={0,1,4,5,7,6,3,2};
    dump(data, "after 4th compare");
    permute(work, data, tab4,8); //4
    dump(work, "after 5th permute");
    compare(work,8);
    INT64T tab5[]={0,4,2,6,1,5,3,7};
    permute(data,work, tab5,8); //5
    compare(data,8);
    INT64T tab6[]={0,2,4,6,1,3,5,7};
    permute(work,data, tab6,8); //6
    copy(data, work, 8);
}

static M256DI store_mask[9][2];
static M256DI load_maxval;
void initialize_lsmask()
{
    int64_t maxval[4];
    int ii;
    for(ii=0;ii<4; ii++){
	maxval[ii]= INT64_MAX;
    }
    load_maxval.d = _mm256_loadu_pd((double*)maxval);
    for(int i=0;i<9;i++){
	uint64_t mask[4];
	for(ii=0;ii<4; ii++){ // a part
	    if( ii<i) {
		mask[ii]= UINT64_MAX;
	    }else{
		mask[ii]= 0;
	    }
	}
	store_mask[i][0].d =_mm256_loadu_pd((double*)mask);
	for(ii=0;ii<4; ii++){ // b part
	    if( ii+4<i) {
		mask[ii]= UINT64_MAX;
	    }else{
		mask[ii]= 0;
	    }
	}
	store_mask[i][1].d =_mm256_loadu_pd((double*)mask);
    }
}


#define compareavx2(a,b,c,e)			\
    {						\
	register M256DI maxflg;			\
	maxflg.i = _mm256_cmpgt_epi64(a.i,b.i);	 \
	c.d = _mm256_blendv_pd(a.d, b.d,  maxflg.d);	\
	e.d = _mm256_blendv_pd(b.d, a.d,  maxflg.d);	\
    }
void bitonic8(INT64T * data, int n)
{
    M256DI a0,b0,a1,b1;
    int64_t imax = INT64_MAX;
    a0.d=_mm256_loadu_pd((double*)data);
    b0.d=_mm256_loadu_pd((double*)(data+4));
    a0.d = _mm256_blendv_pd(load_maxval.d, a0.d,  store_mask[n][0].d);
    b0.d = _mm256_blendv_pd(load_maxval.d, b0.d,  store_mask[n][1].d);
    
    //    dumpavx2(a0,b0, "avx2 content");
    compareavx2(a0,b0,a1,b1);
    //    dumpavx2(a1,b1, "after compare avx2");
    a0 = a1;
    b0.d = _mm256_shuffle_pd(b1.d, b1.d, 0b0101);
    //    dumpavx2(a0,b0, "after second permute avx2");
    compareavx2(a0,b0,a1,b1);
    //    dumpavx2(a1,b1, "after 2nd compare avx2");
    M256DI tmp;
    a0.d = _mm256_shuffle_pd(a1.d, b1.d, 0b1010);
    b0.d = _mm256_shuffle_pd(a1.d, b1.d, 0b0101);
    //    dumpavx2(a0,b0, "after 3rd permute avx2");
    
    compareavx2(a0,b0,a1,b1);
    //    dumpavx2(a1,b1, "after 3rd comp avx2");
    a1.i = _mm256_permute4x64_epi64(a1.i, 0b11011000);
    b1.i = _mm256_permute4x64_epi64(b1.i, 0b01100011);
    a0.d = _mm256_blend_pd(a1.d, b1.d, 0b1010);
    b0.d = _mm256_blend_pd(a1.d, b1.d, 0b0101);
    b0.i = _mm256_permute4x64_epi64(b0.i, 0b01101100);
    //    dumpavx2(a0,b0, "after 4th permute avx2");
    compareavx2(a0,b0,a1,b1);
    //    dumpavx2(a1,b1, "after 4th compare avx2");
    a0.d = _mm256_blend_pd(a1.d, b1.d, 0b1100);
    b0.d = _mm256_blend_pd(a1.d, b1.d, 0b0011);
    a0.i = _mm256_permute4x64_epi64(a0.i, 0b10110100);
    b0.i = _mm256_permute4x64_epi64(b0.i, 0b00011110);
    //    dumpavx2(a0,b0, "after 5th permute avx2");
    compareavx2(a0,b0,a1,b1);
    //    dumpavx2(a1,b1, "after 5th compare avx2");
    b1.i = _mm256_permute4x64_epi64(b1.i,0b10110001);
    a0.d = _mm256_blend_pd(a1.d, b1.d, 0b1010);
    b0.d = _mm256_blend_pd(a1.d, b1.d, 0b0101);
    b0.i = _mm256_permute4x64_epi64(b0.i, 0b10110001);
    //    dumpavx2(a0,b0, "after 6th permute avx2");
    compareavx2(a0,b0,a1,b1);
    //    dumpavx2(a1,b1, "after 6th compare avx2");
    
    b1.i = _mm256_permute4x64_epi64(b1.i,0b01001110);
    a0.d = _mm256_blend_pd(a1.d, b1.d, 0b1100);
    b0.d = _mm256_blend_pd(a1.d, b1.d, 0b0011);
    //    dumpavx2(a0,b0, "after 7th brend   avx2");
    b0.i = _mm256_permute4x64_epi64(b0.i, 0b01110010);
    a0.i = _mm256_permute4x64_epi64(a0.i, 0b11011000);
    //    dumpavx2(a0,b0, "after 7th permute avx2");
    _mm256_maskstore_pd((double*)data,store_mask[n][0].i,  a0.d);
    _mm256_maskstore_pd((double*)(data+4),store_mask[n][1].i,  b0.d);
    //    _mm256_storeu_pd((double*)(data+4), b0.d);
    
}
