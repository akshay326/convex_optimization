/*
reset;model ex2_2_1.mod;data ex2_2_1.dat;solve;
*/


param T;
set states;
set actions;

param NA {actions, states} binary ;	# state-action combinations that aint allowed
param r {actions,states};
param p {states, actions, states};

var f{a in actions,s in states} >= 0, <= NA[a,s];
	# a combination thats not allowed, f = 0

# PTM, probab transition matrix, P(f)
var P{s_next in states, s in states} =
	sum{a in actions} p[s,a,s_next]*f[a,s];

# Immediate expected reward vector, R(f)
var R{s in states} = sum{a in actions} r[a,s]*f[a,s];

var v =
	sum{s_next in states, s in states}	(1 + P[s_next,s] + P[s_next,s]*P[s_next,s]);


minimize objective: 0;

subject to sum_f{s in states}:
	sum{a in actions} f[a,s] = 1;
	
subject to sum_p{s in states, a in actions}:
	sum{s_next in states} p[s,a,s_next] = NA[a,s];

subject to time_pass: f[1,'s1'] = 0;
	