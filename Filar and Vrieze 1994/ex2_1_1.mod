/*
reset;model ex2_1_1.mod;data ex2_1_1.dat;solve;
*/


param T;
param beta = 0.8;
set states;
set actions;

param NA {actions, states} binary ;	# state-action combinations that aint allowed
param r {actions,states};
param p {states, actions, states};

param f{a in actions,s in states};

# PTM, probab transition matrix, P(f)
var P{s_next in states, s in states} =
	sum{a in actions} p[s,a,s_next]*f[a,s];

# Immediate expected reward vector, R(f)
var R{s in states} = sum{a in actions} r[a,s]*f[a,s];

	# Donno how to get inverse, so leaving in this form
var v_inv{s_next in states, s in states}= 
	(if s_next = s then 1 else 0) - beta*P[s_next,s];


minimize objective: 0;

subject to sum_f{s in states}:
	sum{a in actions} f[a,s] = 1;
	
subject to sum_p{s in states, a in actions}:
	sum{s_next in states} p[s,a,s_next] = NA[a,s];
	
