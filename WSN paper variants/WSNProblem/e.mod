# reset;model nominal_solution.mod;data stoch.dat;solve;

param W;
param discount;
var Pd;


set actions = {"low","high"};
set states = {0..W};
set players = {1,2};


param cost {w in states} = 10 * (exp(0.05*w)-1);
param constants{actions,players};
param S {actions,players};


var sigma {players, states, actions};
var value {players, states} ;

#IDS is represented by Player2

param XI {i in players, w in states, a1 in actions, a2 in actions} = 
	if i=1 then
		-cost[w] + constants[a1,1]
	else
		 cost[w] + constants[a2,1];

var defend_probab {a1 in actions, a2 in actions} = min( 0.16*Pd*S[a2,2]/(S[a1,1] + 0.000001), 1);

var PI {w in states, w_next in states, a in actions, b in actions} =  
					(if w_next = w-1 then (w*defend_probab[a,b]/W) 
					
					else (if w_next = w then (w*(1-defend_probab[a,b])/W  + (W-w)*defend_probab[a,b]/W)   
					
					else  (if w_next = w+1 then  ((W-w)*(1-defend_probab[a,b])/W)  else 0   )    ))  ;
					

var C { i in players, w in states} = 
	sum{a1 in actions, a2 in actions}  XI[i,w,a1,a2]*sigma[1,w,a1]*sigma[2,w,a2];


var PIw {w in states, w_next in states} = 
	sum{a1 in actions, a2 in actions}  PI[w,w_next,a1,a2]*sigma[1,w,a1]*sigma[2,w,a2];


minimize total_val : 
	sum {i in players, w in states}(
		 C[i,w] + discount*sum {w_next in states} PIw[w, w_next]*value[i,w_next] - value[i,w] 
	);

subject to indiv_value {i in players, w in states,a in actions}:
	value [i,w] <=  C[i,w] + discount*sum {w_next in states} PIw[w, w_next]*value[i,w_next];
		
subject to probab_sum { i in players, w in states } : 
	sum {a in actions} sigma[i,w,a] = 1 ;

subject to nonnegative {i in players, w in states, a in actions} : 
	sigma[i,w,a]>=0;


