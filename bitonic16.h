void dumpavx512(__m512i a,
                __m512i b,
                char *c)
{
    int64_t val[16];
    _mm512_storeu_si512(val, a);
    _mm512_storeu_si512(val+8, b);
    fprintf(stderr, "%s", c);
    for(int i=0;i<16;i++){fprintf(stderr, " %ld", val[i]);}
    fprintf(stderr, "\n");
}    


#define compare_and_swap16(a, b, a1, b1){    a1 = _mm512_min_epi64(a,b);    b1 = _mm512_max_epi64(a,b);}

   static int64_t __attribute__ ((aligned(64)))  indexc[8]={0,8,1,9,2,10,3,11};
   static int64_t __attribute__ ((aligned(64)))  indexd[8]={4,12,5,13,6,14,7,15};
 
void  final_reorder16(__m512i a, __m512i b, int64_t* data, int  n)
{
   int m= 16>>1;
   __m512i c;
   __m512i d;
   c=  _mm512_permutex2var_epi64 (a, _mm512_load_epi64(indexc), b);
   d=  _mm512_permutex2var_epi64 (a, _mm512_load_epi64(indexd), b);
   if (n<= m){
      	__mmask8 maska =(__mmask8)  ((1<<n)-1);
         _mm512_mask_storeu_epi64(data, maska, c);
   }else{
   __mmask8 maskb =(__mmask8)  ((1<<n)-1);
   _mm512_storeu_si512(data, c);
   _mm512_mask_storeu_epi64(data+m, maskb, d);
   }
}


#define  initial_copy16(data, a, b, n){   int64_t intmax = INT64_MAX;   a= _mm512_broadcastq_epi64(*((__m128i*)(&intmax)));   b= _mm512_broadcastq_epi64(*((__m128i*)(&intmax)));   __mmask8 maska =(__mmask8)  ((1<<n)-1);   __mmask8 maskb =(__mmask8)  (((1<<n)-1)>>8);   a=_mm512_mask_loadu_epi64(a,maska, data);   b=_mm512_mask_loadu_epi64(b,maskb, data+8);}

void bitonic16(int64_t *data, int n)
{
    int m = 8;
    __m512i a, b, a1, b1;
    initial_copy16(data, a, b, n);
//    dumpavx512(a, b, "after initial copy");
    compare_and_swap16(a,b,a1,b1);
//    dumpavx512(a1, b1, "after 1st compare");
    static int64_t __attribute__ ((aligned(64)))  index0[8]={0,8,2,10,4,12,6,14};
    a= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index0), b1);
    static int64_t __attribute__ ((aligned(64)))  index1[8]={9,1,11,3,13,5,15,7};
    b= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index1), b1);

    compare_and_swap16(a,b,a1,b1);
    static int64_t __attribute__ ((aligned(64)))  index2[8]={0,9,2,11,4,13,6,15};
    a= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index2), b1);
    static int64_t __attribute__ ((aligned(64)))  index3[8]={1,8,3,10,5,12,7,14};
    b= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index3), b1);

    compare_and_swap16(a,b,a1,b1);
    static int64_t __attribute__ ((aligned(64)))  index4[8]={0,8,1,9,4,12,5,13};
    a= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index4), b1);
    static int64_t __attribute__ ((aligned(64)))  index5[8]={11,3,10,2,15,7,14,6};
    b= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index5), b1);

    compare_and_swap16(a,b,a1,b1);
    static int64_t __attribute__ ((aligned(64)))  index6[8]={0,1,11,10,4,5,15,14};
    a= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index6), b1);
    static int64_t __attribute__ ((aligned(64)))  index7[8]={2,3,9,8,6,7,13,12};
    b= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index7), b1);

    compare_and_swap16(a,b,a1,b1);
    static int64_t __attribute__ ((aligned(64)))  index8[8]={0,8,2,10,4,12,6,14};
    a= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index8), b1);
    static int64_t __attribute__ ((aligned(64)))  index9[8]={1,9,3,11,5,13,7,15};
    b= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index9), b1);

    compare_and_swap16(a,b,a1,b1);
    static int64_t __attribute__ ((aligned(64)))  index10[8]={0,8,1,9,2,10,3,11};
    a= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index10), b1);
    static int64_t __attribute__ ((aligned(64)))  index11[8]={15,7,14,6,13,5,12,4};
    b= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index11), b1);

    compare_and_swap16(a,b,a1,b1);
    static int64_t __attribute__ ((aligned(64)))  index12[8]={0,1,2,3,15,14,13,12};
    a= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index12), b1);
    static int64_t __attribute__ ((aligned(64)))  index13[8]={4,5,6,7,11,10,9,8};
    b= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index13), b1);

    compare_and_swap16(a,b,a1,b1);
    static int64_t __attribute__ ((aligned(64)))  index14[8]={0,1,8,9,4,5,12,13};
    a= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index14), b1);
    static int64_t __attribute__ ((aligned(64)))  index15[8]={2,3,10,11,6,7,14,15};
    b= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index15), b1);

    compare_and_swap16(a,b,a1,b1);
    static int64_t __attribute__ ((aligned(64)))  index16[8]={0,8,2,10,4,12,6,14};
    a= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index16), b1);
    static int64_t __attribute__ ((aligned(64)))  index17[8]={1,9,3,11,5,13,7,15};
    b= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index17), b1);

    compare_and_swap16(a,b,a1,b1);
{
   __m512i c;
   __m512i d;
   c=  _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(indexc), b1);
   d=  _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(indexd), b1);
   if (n<= m){
      	__mmask8 maska =(__mmask8)  ((1<<n)-1);
         _mm512_mask_storeu_epi64(data, maska, c);
   }else{
   __mmask8 maskb =(__mmask8)  (((1<<n)-1)>>m);
   _mm512_storeu_si512(data, c);
   _mm512_mask_storeu_epi64(data+m, maskb, d);
   }
}
}
