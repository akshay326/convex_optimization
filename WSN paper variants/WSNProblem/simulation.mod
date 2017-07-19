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
	Creating PI and its weighted counterpart, PIw
*/
param Pd := Uniform(Pmin,Pmax);
param p {a in 1..L, b in 1..L} = min(0.16*Pd*S[2,b]/(S[1,a] + EPS), 1); # EPS To prevent 1/0

param PI {w in 0..M, w_next in 0..M, a in 1..L, b in 1..L} =  
	if w_next = w-1 then 
		(w*p[a,b]/M) 
	else if w_next = w then
		(w*(1-p[a,b])/M  + (M-w)*p[a,b]/M)   
	else if w_next = w+1 then  
		(M-w)*(1-p[a,b])/M 
	else 
		0;


var PIw {w1 in 0..M, w2 in 0..M} = 
	sum {a in 1..L, b in 1..L} sigma[1,w1,a]*sigma[2,w1,b]*PI[w1,w2,a,b];


/* 
	Creating XI and its weighted counterpart, C[i,w,s(w)]
	
	Note:
	1) In WSN paper, Xi[i, w, s(w)] where s(w) is a vector of len L^N
		So, XI is N*M*L^N
	2) C[i,w,s(w)] following eqn (10) in WSN paper
		
*/

param cost {w in 0..M} = alpha * (exp(beta*w)-1);

param XI {i in 1..N, w in 0..M, a1 in 1..L, a2 in 1..L } =	
	if i=1 then
		-cost[w] + constants[1,a1]		# XI1 part
	else 
		cost[w] + constants[2,a2];		# XI2 part

/*
	Weighted XI shud be matrix of L^N-1 rows n L cols		
	Ref #17, pg 9, top-right, eqn (8) states XI[a(-i), a(i)]
	The a(i) part can take values b/t 1..L [Columns]
	While a(-i) part takes L^(N-1) combinations [Rows]
	
	Here was the confusion!!
	In #17, this matrix is a function of 'i', 'w'
	While in WSN paper, it depends on s(w) as well - which is wrong 

var C {i in 1..N, w in 0..M,
		a in 1..L,   	# If N=3, it had been a1 in 1..L, a2 in 1..L 
		col in 1..L } =	
		if i = 1 then
			XI[1,w,col,a]*sigma[2,w,col]
		else if i = 2 then
			XI[2,w,a,col]*sigma[1,w,col];
			
#		else if i = 3 then
#		XI[i,w,a1,a2,col]...
*/

var C {i in 1..N, w in 0..M, action in 1..L} =	
	if i = 1 then
		XI[i,w,action,1]*sigma[2,w,action] + XI[i,w,action,2]*sigma[2,w,action]
	else if i = 2 then
		XI[i,w,1,action]*sigma[1,w,action] + XI[i,w,2,action]*sigma[1,w,action];



minimize diff_to_ideal_value: 
		sum {i in 1..N, w in 0..M} (
			value[i,w]
			- (sum {l in 1..L} C[i,w,l]*sigma[i,w,l])
			- (delta*sum {w2 in 0..M} (value[i,w2]*PIw[w2, w]))
		);

subject to max_condition {i in 1..N, w in 0..M, a in 1..L}:
		value[i, w] >= C[i,w,a] 
			   + delta*sum {b in 1..L, w2 in 0..M} value[i,w2]*(if i=1  then sigma[2,w,b]*PI[w, w2,a,b] else PI[w, w2,b,a]*sigma[1,w,b]);

subject to probab_sum {i in 1..N, w in 0..M}:
		sum {a in 1..L} (sigma[i,w,a]) = 1;
		
subject to positive_probab {i in 1..N, w in 0..M, a in 1..L} :
		sigma[i,w,a] >= 0;