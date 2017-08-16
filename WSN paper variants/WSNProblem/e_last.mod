# reset;model e.mod;data stoch.dat;solve;

param W;
param discount;
var Pd = 0.8;
fix Pd;


set actions = {"low","high"};
set states = {0..W};
set players = {1,2};


param cost {w in states} = 10 * (exp(0.05*w)-1);
param constants{actions,players};
param S {actions,players};


var sigma {players, states, actions};
var value {players, states} ;

# inrtuder -> 1 and IDS -> 2 

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
					


var C { i in players, w in states, a in actions} = 
	if i = 1 then
		(sum{b in actions} XI[i,w,a,b])*sigma[i,w,a]
	else
		(sum{b in actions} XI[i,w,b,a])*sigma[i,w,a];


var PIw {w in states, w_next in states} = 
	sum { a in actions, b in actions } sigma[1,w,a]*sigma[2,w,b]*PI[w,w_next,a,b] ;


minimize total_val : 
		sum {i in players, w in states}(
			 value[i,w]
			- (sum {a in actions} C[i,w,a]*sigma[i,w,a])
			- (discount*sum {w2 in states} (value[i,w2]*PIw[w2, w])) 
		);

subject to indiv_value {i in players, w in states,a in actions}:
	value[i, w] >= C[i,w,a] 
	+ discount*sum {b in actions, w2 in states} value[i,w2]*(if i=1 then PI[w,w2,a,b]*sigma[2,w,b] else PI[w,w2,b,a]*sigma[1,w,b]);
		
subject to probab_sum { i in players, w in states } : sum {a in actions} sigma[i,w,a] = 1 ;

subject to nonnegative {i in players, w in states, a in actions} : sigma[i,w,a]>=0;


