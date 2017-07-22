 set N := 1..9 ;
 set DATA within {N cross N cross N};
 var x {(i,j,k) in {N cross N cross N}} binary ;
 
 
var result{i in N, j in N} = sum {k in N} k*x[k,i,j];


 minimize nothing : 0;

 subject to Columns {j in N, k in N}:  sum {i in N} x[i,j,k] = 1;
 subject to Rows 	{i in N, k in N}:  sum {j in N} x[i,j,k] = 1;
 subject to filled  {i in N, j in N}:  sum {k in N} x[i,j,k] = 1;
 subject to known {(i,j,k) in DATA }:               x[i,j,k] = 1;
 
 
 subject to Squares {k in N, p in 1..3 , q in 1..3}:
 	sum {j in (3*p -2) ..(3* p)} sum {i in (3*q -2) ..(3* q)} x[i,j,k] = 1;