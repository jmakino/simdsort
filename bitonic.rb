def bitonic(n,m)
  if n==2
    butterfly(2,m)
  else      
    bitonic(n/2,m)+
    invert_and_compare(n,m)+
    butterfly(n/2,m)
  end
end
def invert_and_compare(n,m)
  "invert_and_compare #{n} #{m} \n"
end
def butterfly(n,m)
  s=""
  while (n>1)
    s+= "butterfly  #{n} #{m} \n"
    n/=2
  end
  s
end

def iac_location_to_comparator(i,n,m)
  iblock = i / n
  ilocal = i%n
  if ilocal <n/2
    icomparator = ilocal
    iwire =0
  else
    icomparator = n-1-ilocal
    iwire =1
  end
  icomparator +=  iblock *(n/2)
  [icomparator, iwire]
end
def iac_comparator_to_location(icomp, iwire,n,m)
  iblock = icomp/(n/2)
  icomplocal = icomp%(n/2)
  iblock*n + icomplocal + iwire*(n-1-icomplocal*2)
end

def butterfly_comparator_to_location(icomp, iwire,n,m)
  iblock = icomp/(n/2)
  icomplocal = icomp%(n/2)
  iblock*n + icomplocal + iwire*(n/2)
end

def butterfly_location_to_comparator(i,n,m)
  iblock = i / n
  ilocal = i%n
  iwire = ilocal/(n/2)
  icomparator = ilocal %(n/2) + iblock *(n/2)
  [icomparator, iwire]
end

def location_to_comparator(type, i,n,m)
  if type == "butterfly"
    butterfly_location_to_comparator(i,n,m)
  else
    iac_location_to_comparator(i,n,m)
  end
end

def comparator_to_location(type, icomp, iwire,n,m)
  if type == "butterfly"
    butterfly_comparator_to_location(icomp, iwire,n,m)
  else
    iac_comparator_to_location(icomp, iwire,n,m)
  end
end

def portstring(icomparator, iwire)
  "ab"[iwire]+(icomparator.to_s)
end

def elementstring(icomparator, iwire, addition)
  "ab"[iwire]+addition+"["+(icomparator.to_s)+"]"
end

def portyloc(icomparator, iwire, n)
  n- (iwire*(n/2)+(icomparator))-0.5
end
def print_swap_code(type0, type1, n1, n2, m)
  2.times{|iwire|
    (m/2).times{|i|
      iloc= comparator_to_location(type1, i, iwire, n2,m)
      src =location_to_comparator(type0, iloc, n1, m) 
      print "    ", elementstring(i,iwire,""), "="
      print elementstring(*src,"1"), ";\n"
    }
  }
  print "\n"
end

def print_swap_code_avx512(type0, type1, n1, n2, m, count)
  2.times{|iwire|
    index=[]
    (m/2).times{|i|
      iloc= comparator_to_location(type1, i, iwire, n2,m)
      src =location_to_comparator(type0, iloc, n1, m)
      index.push src[1]*(m/2)+src[0];
    }
    print "    static int64_t __attribute__ ((aligned(64)))  index#{count*2+iwire}[#{m/2}]={"+
          index.join(",") +
          "};\n"
    print "    "+ "ab"[iwire] + "= _mm512_permutex2var_epi64 (a1, _mm512_load_epi64(index#{count*2+iwire}), b1);\n"
  }
  print "\n"
end

def print_connections(type0, type1, n1, n2, m)
  2.times{|iwire|
    (m/2).times{|i|
      iloc= comparator_to_location(type1, i, iwire, n2,m)
      src =location_to_comparator(type0, iloc, n1, m) 
      print  portstring(*src), " - ", portstring(i,iwire), "\n"
    }
  }
  print "\n"
end

def print_connection_as_wires(type0, type1, n1, n2, m, x,dx)
  2.times{|iwire|
    (m/2).times{|i|
      iloc= comparator_to_location(type1, i, iwire, n2,m)
      src =location_to_comparator(type0, iloc, n1, m)
      print "reloc ", x, " ",  portyloc(*src, m),"\n"
      print "draw ", x+dx, " ",  portyloc(i,iwire,m), "\n"
    }
  }
  print "\n"
end
    

def print_comparators(type0,  n,  m, x,dx)
  (m/2).times{|i|
    iloc0= comparator_to_location(type0, i, 0, n,m)
    iloc1= comparator_to_location(type0, i, 1, n,m)
    print "reloc ", x, " ", m-0.5- iloc0,"\n"
    print "draw ", x, " ",  m-0.5-iloc1, "\n"
    x = x + dx.to_f/(m/2)
  }
  print "\n"
end


def print_comparators_for_sorter(s,m)
  print "printer test.eps/vcps\n"
  print "square\n"
  a=s.chomp.split("\n")
  print "limit -0.5 #{a.size} 0 #{m}\n"
  print "box\n"
  (a.size()).times{|i|
    from=a[i].split
    print_comparators(from[0], from[1].to_i, from[2].to_i, i, 0.7)
  }
end
def print_wires_for_sorter(s,m)
  print "printer test.eps/vcps\n"
  print "square\n"
  a=s.chomp.split("\n")
  print "limit -0.5 #{a.size} 0 #{m}\n"
  print "box\n"
  (a.size()-1).times{|i|
    from=a[i].split
    to=a[i+1].split
    #  print_connection(from[0], to[0], from[1].to_i,to[1].to_i,
    #                   from[2].to_i)
  print_connection_as_wires(from[0], to[0], from[1].to_i,to[1].to_i,
                  from[2].to_i, i, 0.7)
}
end  


def print_connection_for_sorter(s,m)
  a=s.chomp.split("\n")
  (a.size()-1).times{|i|
    from=a[i].split
    to=a[i+1].split
    print_connections(from[0], to[0], from[1].to_i,to[1].to_i,
                     from[2].to_i)
}
end  

def generate_generic_sorter(s,m)

  print <<-EOF  
void compare_and_swap#{m}(int64_t * a,
                       int64_t * b,
		       int64_t * a1,
		       int64_t *b1)
{
   for(int i=0;i< #{m/2}; i++){
       a1[i]= a[i]<b[i]? a[i]:b[i];
       b1[i]= a[i]<b[i]? b[i]:a[i];
   }
}
 
void final_reorder#{m}(int64_t * a,
                   int64_t * b,
                   int64_t * data,
                   int n)
{
   for(int i=0;i<n; i++){
       if (i%2){
           data[i]= b[i/2];
       }else{
          data[i]= a[i/2];
       }

   }	
}

void initial_copy#{m}(int64_t * data,
                   int64_t * a,
                   int64_t * b,
                   int n)
                   
{
   int m = #{m/2};
   for(int i=0;i< n; i++){
       if (i<m){
          a[i]= data[i];
       }else{
          b[i-m]=data[i];
       }

   }
  for(int i=n;i< #{m}; i++){
       if (i<m){
          a[i]= INT64_MAX;
       }else{
          b[i-m]=INT64_MAX;
       }
   }

      
}

void bitonic#{m}(int64_t *data, int n)
{
    int m = #{m/2};
    int64_t a[m];
    int64_t b[m];
    int64_t a1[m];
    int64_t b1[m];
    initial_copy#{m}(data, a, b, n);
    compare_and_swap#{m}(a,b,a1,b1);
EOF

  a=s.chomp.split("\n")
  (a.size()-1).times{|i|
    from=a[i].split
    to=a[i+1].split
    print_swap_code(from[0], to[0], from[1].to_i,to[1].to_i,
                    from[2].to_i)
    print "    compare_and_swap#{m}(a,b,a1,b1);\n"
  }
  print <<-EOF
    final_reorder#{m}(a1,b1,data, n);
}
EOF
end

def generate_avx512_sorter(s,m)
  if m!= 16
    STDERR.print "Currently only 16 is supported for avx512\n"
    exit -1
  end
  print <<-EOF
void dumpavx512(__m512i a,
                __m512i b,
                char *c)
{
    int64_t val[16];
    _mm512_storeu_si512(val, a);
    _mm512_storeu_si512(val+8, b);
    fprintf(stderr, "%s", c);
    for(int i=0;i<16;i++){fprintf(stderr, " %ld", val[i]);}
    fprintf(stderr, "\\n");
}    


#define compare_and_swap#{m}(a, b, a1, b1)\
{\
    a1 = _mm512_min_epi64(a,b);\
    b1 = _mm512_max_epi64(a,b);\
}

   static int64_t __attribute__ ((aligned(64)))  indexc[8]={0,8,1,9,2,10,3,11};
   static int64_t __attribute__ ((aligned(64)))  indexd[8]={4,12,5,13,6,14,7,15};
 
void  final_reorder#{m}(__m512i a, __m512i b, int64_t* data, int  n)
{
   int m= #{m}>>1;
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


#define  initial_copy#{m}(data, a, b, n)\
{\
   int64_t intmax = INT64_MAX;\
   a= _mm512_broadcastq_epi64(*((__m128i*)(&intmax)));\
   b= _mm512_broadcastq_epi64(*((__m128i*)(&intmax)));\
   __mmask8 maska =(__mmask8)  ((1<<n)-1);\
   __mmask8 maskb =(__mmask8)  (((1<<n)-1)>>#{m/2});\
   a=_mm512_mask_loadu_epi64(a,maska, data);\
   b=_mm512_mask_loadu_epi64(b,maskb, data+#{m/2});\
}

void bitonic#{m}(int64_t *data, int n)
{
    int m = #{m/2};
    __m512i a, b, a1, b1;
    initial_copy#{m}(data, a, b, n);
//    dumpavx512(a, b, "after initial copy");
    compare_and_swap#{m}(a,b,a1,b1);
//    dumpavx512(a1, b1, "after 1st compare");
EOF

  a=s.chomp.split("\n")
  (a.size()-1).times{|i|
    from=a[i].split
    to=a[i+1].split
    print_swap_code_avx512(from[0], to[0], from[1].to_i,to[1].to_i,
                    from[2].to_i, i)
#    print "dumpavx512(a, b, \"after #{i} permute\");\n"
    print "    compare_and_swap#{m}(a,b,a1,b1);\n"
#    print "dumpavx512(a1, b1, \"after #{i} compare\");\n"
  }
  print <<-EOF
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
EOF
end


target="GENERIC"
if ARGV[0]== "AVX512"
  target = "AVX512"
  m=16
else   
  m=ARGV[0].to_i
end
mode = ARGV[1].to_i
s= bitonic(m,m)
#print_comparators_for_sorter(s,m)
if mode==0
  print_wires_for_sorter(s,m)
elsif mode==1
  print_comparators_for_sorter(s,m)
elsif mode==2
  print_connection_for_sorter(s,m)
elsif mode==3
  if target == "GENERIC"
    generate_generic_sorter(s,m)
  elsif target == "AVX512"
    generate_avx512_sorter(s,m)
  end
    
end
