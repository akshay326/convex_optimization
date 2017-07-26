#	reset;model corollary.mod;option show_status 1;solve;


param M = 30;	# No of states
param N = 2;	# No of players
param L = 2;	# No of actions , assumed equal for all players
param Pmin := 0.8;
param Pmax := 0.8;

# 1 -> Intruder, 2 -> IDS
param S {1..N,1..L};	# Set f strategies for every player
let S[1,1] := 0;
let S[1,2] := 8;
let S[2,1] := 10;
let S[2,2] := 50;

param constants {1..N, 1..L};
let constants[1,1] := 0;
let constants[1,2] := 8;
let constants[2,1] := 2;
let constants[2,2] := 10;

param delta = 0.9;
param alpha = 10;
param beta = 0.05;
param EPS = 0.00000001;



/*
	Data Part Ends Here
	Model Starts
*/

var sigma {1..N, 0..M, 1..L};
var value {1..N, 0..M};


/*
	Creating PI and XI
*/
param Pd := Uniform(Pmin,Pmax);
param p {a1 in 1..L, a2 in 1..L} = min(0.16*Pd*S[2,a2]/(S[1,a1] + EPS), 1); # EPS To prevent 1/0

param PI {w in 0..M, w_next in 0..M, a1 in 1..L, a2 in 1..L} =  
	if w_next = w-1 then 
		(w*p[a1,a2]/M) 
	else if w_next = w then
		(w*(1-p[a1,a2])/M  + (M-w)*p[a1,a2]/M)   
	else if w_next = w+1 then  
		(M-w)*(1-p[a1,a2])/M 
	else 
		0;

param cost {w in 0..M} = alpha * (exp(beta*w)-1);

param XI {i in 1..N, w in 0..M, a1 in 1..L, a2 in 1..L } =	
	if i=1 then
		-cost[w] + constants[1,a1]		# XI1 part
	else 
		cost[w] + constants[2,a2];		# XI2 part

/*
	C(i,w) is L^N-1 x L
	So C is N*M*( L*L )
*/
var C {i in 1..N, w in 0..M, a1 in 1..L, a2 in 1..L} = 
	if i = 1 then
		XI[1,w,a1,a2]*sigma[2,w,a2]
	else
		XI[2,w,a1,a2]*sigma[1,w,a1];

/*
	Took A(w) as 2*M*L^N x M*L^N
	So A is M*( 2*....)	
*/
param A{w in 0..M, type in 1..2, 
		w1 in 0..M, a11 in 1..L, a12 in 1..L,
		w2 in 0..M, a21 in 1..L, a22 in 1..L}=
	if a11 = a21 and a12 = a22 and w1=w2 then 
		if type = 1 then  
			1
		else
			-1
	else 0;

/*
	Taking b(w) vector as 2*M*L^N
	b should be M*( 2*....)	
*/
param PMIN {a1 in 1..L, a2 in 1..L} = min(0.16*Pmin*S[2,a2]/(S[1,a1] + EPS), 1);
param PMAX {a1 in 1..L, a2 in 1..L} = min(0.16*Pmax*S[2,a2]/(S[1,a1] + EPS), 1);

param b{w in 0..M, type in 1..2, 
		w_next in 0..M, s1 in 1..L, s2 in 1..L}= 
		if type = 1 then  
			if w_next = w-1 then
				(w*PMIN[s1,s2])/M
			else if w_next = w then
				if 2*w <= M then
					(w/M) + PMIN[s1,s2]*(1-(2*w/M))
				else
					(w/M) + PMAX[s1,s2]*(1-(2*w/M))
			else if w_next = w+1 then
				 (1-(w/M))*(1-PMAX[s1,s2])
			else
				0
		else
			if w_next = w-1 then
				-(w*PMAX[s1,s2])/M
			else if w_next = w then
				if 2*w <= M then
					-(w/M) - PMAX[s1,s2]*(1-(2*w/M))
				else
					-(w/M) - PMIN[s1,s2]*(1-(2*w/M))
			else if w_next = w+1 then
				 -(1-(w/M))*(1-PMIN[s1,s2]);

/*
	D(w) is L^N x M*L^N
	So D shud be M*L^N x M*L^N such that
*/
param D{w in 0..M, a1 in 1..L, a2 in 1..L, 
		w_next in 0..M, s1 in 1..L, s2 in 1..L} = 
		if a1 = s1 and a2 = s2 then 
			1
		else 
			0;

/*
	z(i,w) is vector of len M*L^N
	so z shud be N*M*M*L^N
	
	Z(i,w) is matrix - M*L^N x L
	So Z shud be N*M*( M*L^N x L )
	
	No need to define it
	Just add as constrains at last
*/
var z{i in 1..N, w in 0..M, w_next in 0..M, s1 in 1..L, s2 in 1..L}=
	value[i,w_next]*sigma[1,w,s1]*sigma[2,w,s2];
	
var Z{i in 1..N, w in 0..M, w_next in 0..M, s1 in 1..L, s2 in 1..L,l in 1..L};

/*
	q(i,w) is 2*M*L^N  So q is N*M*( 2*M*L^N )
	r(i,w) is L^N 	   So r is N*M*( L^N )
	t(i,w) is M*L^N    So t is N*M*( M*L^N )
	
	t(i, sigma(w)) is vector of len M
	=[sum{s1 in 1..L, s2 in 1..L} sigma[i,w,s1]*sigma[i,w,s2]*t[i,w,w_next,s1,s2]]
		w_next in 0..M
	No need to define these
*/
var q{1..N, 0..M, 1..2, 0..M, 1..L, 1..L};
var r{1..N, 0..M, 1..L, 1..L};
var t{1..N, 0..M, 0..M, 1..L, 1..L};

/*
	# CHECK 2: A(w)*p(w) >= b(w)
	
	for{w in 0..M, type in 1..2, w1 in 0..M, a1 in 1..L, a2 in 1..L}{
		print((sum {w_next in 0..M, s1 in 1..L, s2 in 1..L} A[w,type,w1,a1,a2,w_next,s1,s2]*PI[w_next,w,s1,s2])-b[w,type,w1,a1,a2]);
	}
*/

/*
	# CHECK 1: D(w)*p(w) = 1
	
	for{w in 0..M, a1 in 1..L, a2 in 1..L}{
		print(sum {w_next in 0..M, s1 in 1..L, s2 in 1..L} (D[w,a1,a2,w_next,s1,s2]*PI[w_next,w,s1,s2]));
	}
*/

#/*
minimize difference: 
	sum {i in 1..N, w in 0..M} (
		value[i,w]
		- (sum {s1 in 1..L, s2 in 1..L} C[i,w,s1,s2]*(if i =1 then sigma[1,w,s1] else sigma[2,w,s2]) )
		- (delta*sum {w2 in 0..M} (value[i,w2]*(sum{s1 in 1..L, s2 in 1..L} sigma[1,w,s1]*sigma[2,w,s2]*t[i,w,w2,s1,s2])))
	);


subject to in1 {i in 1..N, w in 0..M, l_dash in 1..L}:
	  (delta*sum {w_next in 0..M, s1 in 1..L, s2 in 1..L} Z[i,w,w_next,s1,s2,l_dash]*t[i,w,w_next,s1,s2])
	+ (sum{s2 in 1..L} ( if i=1 then C[i,w,l_dash,s2] else C[i,w,s2,l_dash] ) )   
	>= value[i, w];
		
		
subject to in2 {i in 1..N, w in 0..M}:
	value[i, w] >= 
	  (sum {s1 in 1..L, s2 in 1..L} C[i,w,s1,s2]*(if i =1 then sigma[1,w,s1] else sigma[2,w,s2]))
	+ (sum{type in 1..2, w_next in 0..M, s1 in 1..L, s2 in 1..L} b[w,type,w_next,s1,s2]*q[i,w,type,w_next,s1,s2])
	+ (sum{s1 in 1..L, s2 in 1..L} r[i,w,s1,s2]);
	
	
subject to in3 {i in 1..N, w in 0..M, w1 in 0..M, a1 in 1..L, a2 in 1..L}:
	(sum{type in 1..2, w2 in 0..M, s1 in 1..L, s2 in 1..L} q[i,w,type,w2,s1,s2]*A[w,type,w2,s1,s2,w1,a1,a2])
   +(sum{s1 in 1..L, s2 in 1..L} r[i,w,s1,s2]*D[w,s1,s2,w1,a1,a2]) 
	>= delta*z[i,w,w1,a1,a2];		# Coz Z*sigma = z 
	
	
subject to in4 {i in 1..N, w in 0..M, type in 1..L, w1 in 0..M, a1 in 1..L, a2 in 1..L} :
	q[i,w,type,w1,a1,a2] <= 0;


subject to in5 {i in 1..N, w in 0..M, a in 1..L} :
	sigma[i,w,a] >= 0;	


subject to in6 {i in 1..N, w in 0..M, type in 1..2, w1 in 0..M, a1 in 1..L, a2 in 1..L} :
	(sum{w2 in 0..M, s1 in 1..L, s2 in 1..L} 
		A[w,type,w1,a1,a2,w2,s1,s2]*t[i,w,w2,s1,s2]) >= b[w,type,w1,a1,a2];	
	

subject to e1 {i in 1..N, w in 0..M}:
	sum {a in 1..L} (sigma[i,w,a]) = 1;

	
subject to e2 {i in 1..N, w in 0..M, w_next in 0..M, s1 in 1..L, s2 in 1..L} :
	sum {l in 1..L} (Z[i,w,w_next,s1,s2,l]*sigma[i,w,l]) = z[i,w,w_next,s1,s2]; 
	
subject to e3 {i in 1..N, w in 0..M, a1 in 1..L, a2 in 1..L} :
	(sum {w_next in 0..M, s1 in 1..L, s2 in 1..L} D[w,a1,a2,w_next,s1,s2]*t[i,w,w_next,s1,s2]) = 1; 
#*/
