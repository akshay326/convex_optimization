param M = 30;	# No of states
param N = 2;	# No of players
param L = 2;	# No of actions , assumed equal for all players

# 1 -> Intruder, 2 -> IDS
param S {1..N,1..L};	# Set f strategies for every player
let S[1,1] := 0;
let S[1,2] := 8;
let S[2,1] := 10;
let S[2,2] := 50;

param C {1..N, 1..L};
let C[1,1] := 0;
let C[1,2] := 8;
let C[2,1] := 2;
let C[2,2] := 10;

param delta = 0.9;
param alpha = 10;
param beta = 0.05;

/*
	Data Part Ends Here
	Model Part Starts
*/

var sigma {1..N, 1..M, 1..L};
var value {1..N, 1..M};

param PI {1..M,1..M, 1..L, 1..L};	# Note if N =3, then 1..L thrice at end
param p;

/*
	Creating PI and its weighted counterpart, PIw
*/
param Pd := Uniform(0.8,0.8); # TODO : Change it later

for {s1 in 1..L, s2 in 1..L}{
	let p := (if S[1,s1] =0 then 1 else min( 0.16*Pd*S[2,s2]/S[1,s1], 1));
	
	for {i in 1..M, j in 1..M}{
		if j == i - 1 then
        	let PI[i,j,s1,s2] := i / M * p;
        else if j == i then
        	let PI[i,j,s1,s2] := i / M * (1 - p) + (1 - i / M) * p;
    	else if j == i + 1 then 
        	let PI[i,j,s1,s2] := (1 - i / M) * p;
        else
        	let PI[i,j,s1,s2] := 0;
	}
}

var PIw {w1 in 1..M, w2 in 1..M} = 
	sum {a in 1..L, b in 1..L} 
		(sigma[1,w1,a]*sigma[2,w1,b]*PI[w1,w2,a,b]);


/* 
	Creating XI and its weighted counterpart, XIw
	
	Note:
	1) In WSN paper, Xi[i, w, s(w)] where s(w) is a vector of len L^N
	2) Creating XIw or C[i,w,s(w)] as eqn (1.2) in ref #36
	   Not following eqn (10) in WSN paper
		
*/

param XI {i in 1..N, w in 1..M, a in 1..L} =
	((-1)^i)*(alpha*(exp(w*beta) - 1)) + C[i,a];

var XIw {i in 1..N, w in 1..M} =	
	sum {a in 1..L} (sigma[i,w,a]*XI[i,w,a]);

/*

	Optimization Part
	
*/

minimize diff_to_ideal_value: 
		sum {i in 1..N, w in 1..M} 
			(value[i,w] - XIw[i,w] - delta*sum {w2 in 1..M} (value[i,w2]*PIw[w2, w]));
			
subject to max_condition {i in 1..N, w in 1..M, a in 1..L}:
		value[i, w] >= XI[i,w,a] + delta*sum {w2 in 1..M ,b in 1..L} 
								(value[i,w2]*PI[w, w2,a,b]*(if i=1 then sigma[2,w,b] else sigma[1,w,b]));

subject to probab_sum {i in 1..N, w in 1..M}:
		sum {a in 1..L} (sigma[i,w,a]) == 1;
		
subject to positive_probab {i in 1..N, w in 1..M, a in 1..L} :
		sigma[i,w,a] >= 0;