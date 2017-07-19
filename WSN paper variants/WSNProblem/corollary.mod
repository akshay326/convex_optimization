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

var C {i in 1..N, w in 0..M, a1 in 1..L, a2 in 1..L} = 
	if i = 1 then
		XI[1,w,a1,a2]*sigma[2,w,a1]
	else
		XI[2,w,a1,a2]*sigma[1,w,a2];
		

var r{i in 1..N, w in 0..M,l in 1..L^N};

var q{i in 1..N, w in 0..M,l in 1..2*M*L^N};

/*
	Dumb notation :)
	But will work	
*/
param A{row in 1..2*M*L^N, col in 1..M*L^N}=
	if 2*col = row then -1
	else if (2*col-1) = row then 1
	else 0;

/*
	D(w) is L^N x M*L^N
	So D shud be M*L^N x M*L^N such that
	for{w in 0..M, a1 in 1..L, a2 in 1..L}{
		sum {w_next in 0..M, s1 in 1..L, s2 in 1..L} 
			(D[w,a1,a2,w_next,s1,s2]*PI[w_next,w,s1,s2]) = 1
	}
	
	I thought of removing the a1 and a2. But dimensions of 'T' wont satisfy  
*/
param D{w in 0..M, a1 in 1..L, a2 in 1..L, 
		w_next in 0..M, s1 in 1..L, s2 in 1..L} = 
		if a1 = s1 and a2 = s2 then 
			1
		else 
			0;
			


# Plz add the code for b(w)

/*
	z(i,w) is vector of len M*L^N
	so z shud be N*M*M*L^N
*/
var z{i in 1..N, w in 0..M, w_next in 0..M, s1 in 1..L, s2 in 1..L}=
	value[i,w_next]*sigma[1,w,s1]*sigma[2,w,s2];
	
/*
	Z(i,w) is matrix - M*L^N x L
	So Z shud be N*M*( M*L^N x L )
	
	No need to define it
	Just add as constrains at last
*/
var Z{i in 1..N, w in 0..M, w_next in 0..M, s1 in 1..L, s2 in 1..L,l in 1..L};

/*
	q(i,w) is 2*M*L^N
	So q is N*M*( 2*M*L^N )
	
	r(i,w) is L^N
	So r is N*M*( L^N )
*/
var q{1..N,0..M, 2*M*L^N};
var r{1..N,0..M, L^N};

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
		
subject to zConstrain {i in 1..N, w in 0..M, w_next in 0..M, s1 in 1..L, s2 in 1..L} :
		Z[i,w,w_next,s1,s2, 1]*sigma[i,w,1] + Z[i,w,w_next,s1,s2, 2]*sigma[i,w,2] = z[i,w,w_next,s1,s2]; 