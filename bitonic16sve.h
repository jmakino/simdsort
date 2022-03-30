void dumpsve2(svint64_t a,
                svint64_t b,
                char *c)
{
    int64_t val[16];
    svst1_s64(svptrue_b64(),val, a);
    svst1_s64(svptrue_b64(),val+8, b);
    fprintf(stderr, "%s", c);
    for(int i=0;i<16;i++){fprintf(stderr, " %ld", val[i]);}
    fprintf(stderr, "\n");
}    


#define compare_and_swap16(a, b, a1, b1){    a1 = svmin_s64_m(ptrue,a,b);    b1 = svmax_s64_m(ptrue,a,b);}

   static int64_t __attribute__ ((aligned(64)))  indexc[8]={0,8,1,9,2,10,3,11};
   static int64_t __attribute__ ((aligned(64)))  indexd[8]={4,12,5,13,6,14,7,15};
 


#define  initial_copy16(data, a, b, n){   int64_t intmax = INT64_MAX;   svint64_t svintmax= svdup_s64(intmax);   svbool_t maska= svwhilelt_b64(0,n);   svbool_t maskb= svwhilelt_b64(8,n);   a= svsel_s64(maska,svld1_s64(ptrue, data),svintmax);   b= svsel_s64(maskb,svld1_s64(ptrue, data+8),svintmax);}

void bitonic16(int64_t *data, int n)
{
    int m = 8;
    svuint64_t svzero= svdup_u64(0);
    svbool_t sel;
    svbool_t ptrue =svptrue_b64();
    svint64_t a, b, a1, b1;
    initial_copy16(data, a, b, n);
//    dumpsve2(a, b, "after initial copy");
    compare_and_swap16(a,b,a1,b1);
//    dumpsve2(a1, b1, "after 1st compare");
    static uint64_t __attribute__ ((aligned(64)))  index0[8]={0,1,2,3,4,5,6,7};
    static uint64_t __attribute__ ((aligned(64)))  index1[8]={1,0,3,2,5,4,7,6};
    static uint64_t __attribute__ ((aligned(64)))  index2[8]={1,0,1,0,1,0,1,0};
    static uint64_t __attribute__ ((aligned(64)))  index3[8]={0,1,2,3,4,5,6,7};
    a=a1;
    b=svtbl_s64(b1, svld1_u64(ptrue,index1));
    sel = svcmpne_u64(ptrue,svzero,svld1_u64(ptrue,index2));
    a1=svsel_s64(sel, a,b);
    b1=svsel_s64(sel, b,a);
    b=svtbl_s64(b1, svld1_u64(ptrue,index3));
    a=a1;
    compare_and_swap16(a,b,a1,b1);
    static uint64_t __attribute__ ((aligned(64)))  index6[8]={0,1,2,3,4,5,6,7};
    static uint64_t __attribute__ ((aligned(64)))  index7[8]={0,1,2,3,4,5,6,7};
    static uint64_t __attribute__ ((aligned(64)))  index8[8]={1,0,1,0,1,0,1,0};
    static uint64_t __attribute__ ((aligned(64)))  index9[8]={1,0,3,2,5,4,7,6};
    a=a1;
    b=b1;
    sel = svcmpne_u64(ptrue,svzero,svld1_u64(ptrue,index8));
    a1=svsel_s64(sel, a,b);
    b1=svsel_s64(sel, b,a);
    b=svtbl_s64(b1, svld1_u64(ptrue,index9));
    a=a1;
    compare_and_swap16(a,b,a1,b1);
    static uint64_t __attribute__ ((aligned(64)))  index12[8]={0,3,1,2,4,7,5,6};
    static uint64_t __attribute__ ((aligned(64)))  index13[8]={3,0,2,1,7,4,6,5};
    static uint64_t __attribute__ ((aligned(64)))  index14[8]={1,0,1,0,1,0,1,0};
    static uint64_t __attribute__ ((aligned(64)))  index15[8]={0,1,2,3,4,5,6,7};
    a=svtbl_s64(a1, svld1_u64(ptrue,index12));
    b=svtbl_s64(b1, svld1_u64(ptrue,index13));
    sel = svcmpne_u64(ptrue,svzero,svld1_u64(ptrue,index14));
    a1=svsel_s64(sel, a,b);
    b1=svsel_s64(sel, b,a);
    b=svtbl_s64(b1, svld1_u64(ptrue,index15));
    a=a1;
    compare_and_swap16(a,b,a1,b1);
    static uint64_t __attribute__ ((aligned(64)))  index18[8]={0,1,2,3,4,5,6,7};
    static uint64_t __attribute__ ((aligned(64)))  index19[8]={1,0,3,2,5,4,7,6};
    static uint64_t __attribute__ ((aligned(64)))  index20[8]={1,1,0,0,1,1,0,0};
    static uint64_t __attribute__ ((aligned(64)))  index21[8]={2,3,0,1,6,7,4,5};
    a=a1;
    b=svtbl_s64(b1, svld1_u64(ptrue,index19));
    sel = svcmpne_u64(ptrue,svzero,svld1_u64(ptrue,index20));
    a1=svsel_s64(sel, a,b);
    b1=svsel_s64(sel, b,a);
    b=svtbl_s64(b1, svld1_u64(ptrue,index21));
    a=a1;
    compare_and_swap16(a,b,a1,b1);
    static uint64_t __attribute__ ((aligned(64)))  index24[8]={0,1,2,3,4,5,6,7};
    static uint64_t __attribute__ ((aligned(64)))  index25[8]={1,0,3,2,5,4,7,6};
    static uint64_t __attribute__ ((aligned(64)))  index26[8]={1,0,1,0,1,0,1,0};
    static uint64_t __attribute__ ((aligned(64)))  index27[8]={1,0,3,2,5,4,7,6};
    a=a1;
    b=svtbl_s64(b1, svld1_u64(ptrue,index25));
    sel = svcmpne_u64(ptrue,svzero,svld1_u64(ptrue,index26));
    a1=svsel_s64(sel, a,b);
    b1=svsel_s64(sel, b,a);
    b=svtbl_s64(b1, svld1_u64(ptrue,index27));
    a=a1;
    compare_and_swap16(a,b,a1,b1);
    static uint64_t __attribute__ ((aligned(64)))  index30[8]={0,7,1,6,2,5,3,4};
    static uint64_t __attribute__ ((aligned(64)))  index31[8]={7,0,6,1,5,2,4,3};
    static uint64_t __attribute__ ((aligned(64)))  index32[8]={1,0,1,0,1,0,1,0};
    static uint64_t __attribute__ ((aligned(64)))  index33[8]={0,1,2,3,4,5,6,7};
    a=svtbl_s64(a1, svld1_u64(ptrue,index30));
    b=svtbl_s64(b1, svld1_u64(ptrue,index31));
    sel = svcmpne_u64(ptrue,svzero,svld1_u64(ptrue,index32));
    a1=svsel_s64(sel, a,b);
    b1=svsel_s64(sel, b,a);
    b=svtbl_s64(b1, svld1_u64(ptrue,index33));
    a=a1;
    compare_and_swap16(a,b,a1,b1);
    static uint64_t __attribute__ ((aligned(64)))  index36[8]={0,1,2,3,4,5,6,7};
    static uint64_t __attribute__ ((aligned(64)))  index37[8]={3,2,1,0,7,6,5,4};
    static uint64_t __attribute__ ((aligned(64)))  index38[8]={1,1,1,1,0,0,0,0};
    static uint64_t __attribute__ ((aligned(64)))  index39[8]={4,5,6,7,0,1,2,3};
    a=a1;
    b=svtbl_s64(b1, svld1_u64(ptrue,index37));
    sel = svcmpne_u64(ptrue,svzero,svld1_u64(ptrue,index38));
    a1=svsel_s64(sel, a,b);
    b1=svsel_s64(sel, b,a);
    b=svtbl_s64(b1, svld1_u64(ptrue,index39));
    a=a1;
    compare_and_swap16(a,b,a1,b1);
    static uint64_t __attribute__ ((aligned(64)))  index42[8]={0,1,2,3,4,5,6,7};
    static uint64_t __attribute__ ((aligned(64)))  index43[8]={2,3,0,1,6,7,4,5};
    static uint64_t __attribute__ ((aligned(64)))  index44[8]={1,1,0,0,1,1,0,0};
    static uint64_t __attribute__ ((aligned(64)))  index45[8]={2,3,0,1,6,7,4,5};
    a=a1;
    b=svtbl_s64(b1, svld1_u64(ptrue,index43));
    sel = svcmpne_u64(ptrue,svzero,svld1_u64(ptrue,index44));
    a1=svsel_s64(sel, a,b);
    b1=svsel_s64(sel, b,a);
    b=svtbl_s64(b1, svld1_u64(ptrue,index45));
    a=a1;
    compare_and_swap16(a,b,a1,b1);
    static uint64_t __attribute__ ((aligned(64)))  index48[8]={0,1,2,3,4,5,6,7};
    static uint64_t __attribute__ ((aligned(64)))  index49[8]={1,0,3,2,5,4,7,6};
    static uint64_t __attribute__ ((aligned(64)))  index50[8]={1,0,1,0,1,0,1,0};
    static uint64_t __attribute__ ((aligned(64)))  index51[8]={1,0,3,2,5,4,7,6};
    a=a1;
    b=svtbl_s64(b1, svld1_u64(ptrue,index49));
    sel = svcmpne_u64(ptrue,svzero,svld1_u64(ptrue,index50));
    a1=svsel_s64(sel, a,b);
    b1=svsel_s64(sel, b,a);
    b=svtbl_s64(b1, svld1_u64(ptrue,index51));
    a=a1;
    compare_and_swap16(a,b,a1,b1);
   svint64_t c;
   svint64_t d;
   c=  svzip1_s64(a1, b1);
   d=  svzip2_s64(a1, b1);
//   dumpsve2(c,d,"after sort and interleave");
   svbool_t maska= svwhilelt_b64(0,n);
   svbool_t maskb= svwhilelt_b64(8,n);
   svst1_s64(maska, data, c);
   svst1_s64(maskb, data+m, d);
}
