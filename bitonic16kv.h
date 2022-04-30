
#undef  compare_and_swap16
#define compare_and_swap16(au, al, ai, bu, bl, bi, au1, al1, ai1,  bu1, bl1, bi1){    __mmask8 gtmask  = _mm512_cmpgt_epu64_mask(au,bu)|              (_mm512_cmpeq_epu64_mask(au,bu)&_mm512_cmpgt_epu64_mask(al,bl))|              (_mm512_cmpeq_epu64_mask(au,bu)&_mm512_cmpeq_epu64_mask(al,bl)&_mm512_cmpgt_epu64_mask(ai,bi));    bu1 = _mm512_mask_blend_epi64(gtmask, bu, au);    bl1 = _mm512_mask_blend_epi64(gtmask, bl, al);    bi1 = _mm512_mask_blend_epi64(gtmask, bi, ai);    au1 = _mm512_mask_blend_epi64(gtmask, au, bu);    al1 = _mm512_mask_blend_epi64(gtmask, al, bl);    ai1 = _mm512_mask_blend_epi64(gtmask, ai, bi);}

//   static int64_t __attribute__ ((aligned(64)))  indexc[8]={0,8,1,9,2,10,3,11};
//   static int64_t __attribute__ ((aligned(64)))  indexd[8]={4,12,5,13,6,14,7,15};
 
void  final_reorder16(__m512i au, __m512i al, __m512i ai, 
                        __m512i bu, __m512i bl, __m512i bi, 
                        int64_t* datau,int64_t* datal,int64_t* datai,
                        int  n)
{
   int m= 16>>1;
   __m512i cu, cl, ci;
   __m512i du, dl, di;
   cu=  _mm512_permutex2var_epi64 (au, _mm512_load_epi64(indexc), bu);
   cl=  _mm512_permutex2var_epi64 (al, _mm512_load_epi64(indexc), bl);
   ci=  _mm512_permutex2var_epi64 (ai, _mm512_load_epi64(indexc), bi);
   du=  _mm512_permutex2var_epi64 (au, _mm512_load_epi64(indexd), bu);
   dl=  _mm512_permutex2var_epi64 (al, _mm512_load_epi64(indexd), bl);
   di=  _mm512_permutex2var_epi64 (ai, _mm512_load_epi64(indexd), bi);

   if (n<= m){
      	__mmask8 maska =(__mmask8)  ((1<<n)-1);
         _mm512_mask_storeu_epi64(datau, maska, cu);
         _mm512_mask_storeu_epi64(datal, maska, cl);
         _mm512_mask_storeu_epi64(datai, maska, ci);
   }else{
   __mmask8 maskb =(__mmask8)  ((1<<n)-1);
   _mm512_storeu_si512(datau, cu);
   _mm512_storeu_si512(datal, cl);
   _mm512_storeu_si512(datai, ci);
   _mm512_mask_storeu_epi64(datau+m, maskb, du);
   _mm512_mask_storeu_epi64(datal+m, maskb, dl);
   _mm512_mask_storeu_epi64(datai+m, maskb, di);
   }
}


#undef  initial_copy16
#define  initial_copy16(datau, datal, datai,  au, al, ai, bu, bl, bi, n){   M128IA intmax;intmax.i[0]= INT64_MAX;   __m512i datamax= _mm512_broadcastq_epi64(intmax.m);   __mmask8 maska =(__mmask8)  ((1<<n)-1);   __mmask8 maskb =(__mmask8)  (((1<<n)-1)>>8);   au=_mm512_mask_loadu_epi64(datamax, maska, datau);   al=_mm512_mask_loadu_epi64(datamax, maska, datal);   ai=_mm512_mask_loadu_epi64(datamax,maska, datai);   bu=_mm512_mask_loadu_epi64(datamax,maskb, datau+8);   bl=_mm512_mask_loadu_epi64(datamax,maskb, datal+8);   bi=_mm512_mask_loadu_epi64(datamax,maskb, datai+8);}

void bitonic16(uint64_t *datau, uint64_t *datal, uint64_t *datai, int n)
{
    int m = 8;
    __m512i au, al, ai, bu, bl, bi, au1, al1, ai1,  bu1, bl1, bi1;
    initial_copy16(datau, datal, datai, au, al, ai, bu, bl, bi, n);
//    dumpavx512(au, bu, "u after initial copy");
//    dumpavx512(al, bl, "l after initial copy");
//    dumpavx512(ai, bi, "i after initial copy");
    compare_and_swap16(au, al, ai, bu, bl, bi, au1, al1, ai1,  bu1, bl1, bi1);
//    dumpavx512(a1, b1, "after 1st compare");
    __m512i tab;
    static int64_t __attribute__ ((aligned(64)))  index0[8]={0,8,2,10,4,12,6,14};
    tab = _mm512_load_epi64(index0);
    au= _mm512_permutex2var_epi64 (au1, tab, bu1);
    al= _mm512_permutex2var_epi64 (al1, tab, bl1);
    ai= _mm512_permutex2var_epi64 (ai1, tab, bi1);
    static int64_t __attribute__ ((aligned(64)))  index1[8]={9,1,11,3,13,5,15,7};
    tab = _mm512_load_epi64(index1);
    bu= _mm512_permutex2var_epi64 (au1, tab, bu1);
    bl= _mm512_permutex2var_epi64 (al1, tab, bl1);
    bi= _mm512_permutex2var_epi64 (ai1, tab, bi1);

    compare_and_swap16(au, al, ai, bu, bl, bi, au1, al1, ai1,  bu1, bl1, bi1);
    static int64_t __attribute__ ((aligned(64)))  index2[8]={0,9,2,11,4,13,6,15};
    tab = _mm512_load_epi64(index2);
    au= _mm512_permutex2var_epi64 (au1, tab, bu1);
    al= _mm512_permutex2var_epi64 (al1, tab, bl1);
    ai= _mm512_permutex2var_epi64 (ai1, tab, bi1);
    static int64_t __attribute__ ((aligned(64)))  index3[8]={1,8,3,10,5,12,7,14};
    tab = _mm512_load_epi64(index3);
    bu= _mm512_permutex2var_epi64 (au1, tab, bu1);
    bl= _mm512_permutex2var_epi64 (al1, tab, bl1);
    bi= _mm512_permutex2var_epi64 (ai1, tab, bi1);

    compare_and_swap16(au, al, ai, bu, bl, bi, au1, al1, ai1,  bu1, bl1, bi1);
    static int64_t __attribute__ ((aligned(64)))  index4[8]={0,8,1,9,4,12,5,13};
    tab = _mm512_load_epi64(index4);
    au= _mm512_permutex2var_epi64 (au1, tab, bu1);
    al= _mm512_permutex2var_epi64 (al1, tab, bl1);
    ai= _mm512_permutex2var_epi64 (ai1, tab, bi1);
    static int64_t __attribute__ ((aligned(64)))  index5[8]={11,3,10,2,15,7,14,6};
    tab = _mm512_load_epi64(index5);
    bu= _mm512_permutex2var_epi64 (au1, tab, bu1);
    bl= _mm512_permutex2var_epi64 (al1, tab, bl1);
    bi= _mm512_permutex2var_epi64 (ai1, tab, bi1);

    compare_and_swap16(au, al, ai, bu, bl, bi, au1, al1, ai1,  bu1, bl1, bi1);
    static int64_t __attribute__ ((aligned(64)))  index6[8]={0,1,11,10,4,5,15,14};
    tab = _mm512_load_epi64(index6);
    au= _mm512_permutex2var_epi64 (au1, tab, bu1);
    al= _mm512_permutex2var_epi64 (al1, tab, bl1);
    ai= _mm512_permutex2var_epi64 (ai1, tab, bi1);
    static int64_t __attribute__ ((aligned(64)))  index7[8]={2,3,9,8,6,7,13,12};
    tab = _mm512_load_epi64(index7);
    bu= _mm512_permutex2var_epi64 (au1, tab, bu1);
    bl= _mm512_permutex2var_epi64 (al1, tab, bl1);
    bi= _mm512_permutex2var_epi64 (ai1, tab, bi1);

    compare_and_swap16(au, al, ai, bu, bl, bi, au1, al1, ai1,  bu1, bl1, bi1);
    static int64_t __attribute__ ((aligned(64)))  index8[8]={0,8,2,10,4,12,6,14};
    tab = _mm512_load_epi64(index8);
    au= _mm512_permutex2var_epi64 (au1, tab, bu1);
    al= _mm512_permutex2var_epi64 (al1, tab, bl1);
    ai= _mm512_permutex2var_epi64 (ai1, tab, bi1);
    static int64_t __attribute__ ((aligned(64)))  index9[8]={1,9,3,11,5,13,7,15};
    tab = _mm512_load_epi64(index9);
    bu= _mm512_permutex2var_epi64 (au1, tab, bu1);
    bl= _mm512_permutex2var_epi64 (al1, tab, bl1);
    bi= _mm512_permutex2var_epi64 (ai1, tab, bi1);

    compare_and_swap16(au, al, ai, bu, bl, bi, au1, al1, ai1,  bu1, bl1, bi1);
    static int64_t __attribute__ ((aligned(64)))  index10[8]={0,8,1,9,2,10,3,11};
    tab = _mm512_load_epi64(index10);
    au= _mm512_permutex2var_epi64 (au1, tab, bu1);
    al= _mm512_permutex2var_epi64 (al1, tab, bl1);
    ai= _mm512_permutex2var_epi64 (ai1, tab, bi1);
    static int64_t __attribute__ ((aligned(64)))  index11[8]={15,7,14,6,13,5,12,4};
    tab = _mm512_load_epi64(index11);
    bu= _mm512_permutex2var_epi64 (au1, tab, bu1);
    bl= _mm512_permutex2var_epi64 (al1, tab, bl1);
    bi= _mm512_permutex2var_epi64 (ai1, tab, bi1);

    compare_and_swap16(au, al, ai, bu, bl, bi, au1, al1, ai1,  bu1, bl1, bi1);
    static int64_t __attribute__ ((aligned(64)))  index12[8]={0,1,2,3,15,14,13,12};
    tab = _mm512_load_epi64(index12);
    au= _mm512_permutex2var_epi64 (au1, tab, bu1);
    al= _mm512_permutex2var_epi64 (al1, tab, bl1);
    ai= _mm512_permutex2var_epi64 (ai1, tab, bi1);
    static int64_t __attribute__ ((aligned(64)))  index13[8]={4,5,6,7,11,10,9,8};
    tab = _mm512_load_epi64(index13);
    bu= _mm512_permutex2var_epi64 (au1, tab, bu1);
    bl= _mm512_permutex2var_epi64 (al1, tab, bl1);
    bi= _mm512_permutex2var_epi64 (ai1, tab, bi1);

    compare_and_swap16(au, al, ai, bu, bl, bi, au1, al1, ai1,  bu1, bl1, bi1);
    static int64_t __attribute__ ((aligned(64)))  index14[8]={0,1,8,9,4,5,12,13};
    tab = _mm512_load_epi64(index14);
    au= _mm512_permutex2var_epi64 (au1, tab, bu1);
    al= _mm512_permutex2var_epi64 (al1, tab, bl1);
    ai= _mm512_permutex2var_epi64 (ai1, tab, bi1);
    static int64_t __attribute__ ((aligned(64)))  index15[8]={2,3,10,11,6,7,14,15};
    tab = _mm512_load_epi64(index15);
    bu= _mm512_permutex2var_epi64 (au1, tab, bu1);
    bl= _mm512_permutex2var_epi64 (al1, tab, bl1);
    bi= _mm512_permutex2var_epi64 (ai1, tab, bi1);

    compare_and_swap16(au, al, ai, bu, bl, bi, au1, al1, ai1,  bu1, bl1, bi1);
    static int64_t __attribute__ ((aligned(64)))  index16[8]={0,8,2,10,4,12,6,14};
    tab = _mm512_load_epi64(index16);
    au= _mm512_permutex2var_epi64 (au1, tab, bu1);
    al= _mm512_permutex2var_epi64 (al1, tab, bl1);
    ai= _mm512_permutex2var_epi64 (ai1, tab, bi1);
    static int64_t __attribute__ ((aligned(64)))  index17[8]={1,9,3,11,5,13,7,15};
    tab = _mm512_load_epi64(index17);
    bu= _mm512_permutex2var_epi64 (au1, tab, bu1);
    bl= _mm512_permutex2var_epi64 (al1, tab, bl1);
    bi= _mm512_permutex2var_epi64 (ai1, tab, bi1);

    compare_and_swap16(au, al, ai, bu, bl, bi, au1, al1, ai1,  bu1, bl1, bi1);
  {
   __m512i cu, cl, ci;
   __m512i du, dl, di;
   cu=  _mm512_permutex2var_epi64 (au1, _mm512_load_epi64(indexc), bu1);
   du=  _mm512_permutex2var_epi64 (au1, _mm512_load_epi64(indexd), bu1);
   cl=  _mm512_permutex2var_epi64 (al1, _mm512_load_epi64(indexc), bl1);
   dl=  _mm512_permutex2var_epi64 (al1, _mm512_load_epi64(indexd), bl1);
   ci=  _mm512_permutex2var_epi64 (ai1, _mm512_load_epi64(indexc), bi1);
   di=  _mm512_permutex2var_epi64 (ai1, _mm512_load_epi64(indexd), bi1);
   if (n<= m){
      	__mmask8 maska =(__mmask8)  ((1<<n)-1);
                 _mm512_mask_storeu_epi64(datau, maska, cu);
                 _mm512_mask_storeu_epi64(datal, maska, cl);
                 _mm512_mask_storeu_epi64(datai, maska, ci);

   }else{
   __mmask8 maskb =(__mmask8)  (((1<<n)-1)>>m);
   _mm512_storeu_si512(datau, cu);
   _mm512_mask_storeu_epi64(datau+m, maskb, du);
   _mm512_storeu_si512(datal, cl);
   _mm512_mask_storeu_epi64(datal+m, maskb, dl);
   _mm512_storeu_si512(datai, ci);
   _mm512_mask_storeu_epi64(datai+m, maskb, di);
  
    }
  }
}
