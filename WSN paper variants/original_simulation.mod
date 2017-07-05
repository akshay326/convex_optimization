param M = 30;	# No of states
param N = 2;	# No of players
param L = 2;	# No of actions , assumed equal for all players

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

#param PI {1..M,1..M, 1..L};
#param PIw {1..M,1..M};
param PI {1..M,1..M};
param p;

/*
	First Create PI
	Then its weighted counterpart, PIw
*/
param Pd := Uniform(0.8,0.8); # TODO : Change it later

for {s1 in 1..L, s2 in 1..L}{	# Actions for players
	let p := min(0.16*Pd*s2/s1,1);
	
	for {i in 1..M, j in 1..M}{	# Possible States for players
		if j == i - 1 then
        	let PI[i,j] := i / M * p;
        else if j == i then
        	let PI[i,j] := i / M * (1 - p) + (1 - i / M) * p;
    	else if j == i + 1 then 
        	let PI[i,j] := (1 - i / M) * p;
        else
        	let PI[i,j] := 0;
	}
}

/* 
	Weighted PI is ready
	Simlarly create XI first
	Later create its weighted counterpart, XIw
	
	Note that in WSN paper,
		Xi[i, w, s(w)] where s(w) is a vector of [s(i,w)]
*/

param XI {1..N, 1..M, 1..L};
param XIw {1..N, 1..M};
param cost;

for {i in 1..N, w in 1..M, a in 1..L}{	
	let cost :=  alpha*(exp(w*beta) - 1);
	let XI[i,w,a] := ((-1)^i)*cost + C[i,a];
}

for {i in 1..N, w in 1..M}{	
	let XIw[i,w] := sum {a in 1..L} (sigma[i,w,a]*XI[i,w,a]);
}

/*
	XI and XIw (Weighted) created
	Now Optimization Part
*/

minimize diff_to_ideal_value: 
		sum {i in 1..N, w in 1..M} 
			(value[i,w] - XIw[i,w] - delta*sum {w2 in 1..M} (value[i,w2]*PI[w2, w]));
			
subject to max_condition :
		for {i in 1..N, w in 1..M, a in 1..L}{
			value[i, w] >= XI[i,w,a] + delta*sum {w2 in 1..M}(value[i, w2]*PI[w2,w]);
		}

subject to probab_sum :
		for {i in 1..N, w in 1..M}{
			sum {a in 1..L} sigma[i,w,a] = 1;
		}

subject to positive_probab :
		for {i in 1..N, w in 1..M, a in 1..L}{
			sigma[i,w,a] >= 0;
		}