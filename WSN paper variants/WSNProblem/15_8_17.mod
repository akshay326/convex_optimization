# reset;model 15_8_17.mod;data 15_8_17.dat;solve;

param W;
param discount;
param Pd = 0.8;


set actions = {"low","high"};
set states = {0..W};
set players = {1,2};

# intruder -> 1 and IDS -> 2

param cost {w in states} = 10*(exp(0.05*w)-1);
param C1{actions};
param C2{actions};
param S1{actions};
param S2{actions};


var sigma {players, states, actions} >= 0;
var value {players, states} ;

# Making XI a fn of a1,a2 creates confusion
# In papers, XI1 and XI2 are defined separately
param XI1{w in states, a in actions} = -cost[w] + C1[a];
param XI2{w in states, a in actions} = cost[w] + C2[a];

param defend_probab {a1 in actions, a2 in actions} = min( 0.16*Pd*S2[a2]/(S1[a1] + 0.000001), 1);

param PI {w in states, w_next in states, a in actions, b in actions} =  
	if w_next = w-1 then 
		(w*defend_probab[a,b]/W)
	else if w_next = w then 
		(w*(1-defend_probab[a,b])/W  + (W-w)*defend_probab[a,b]/W)
	else if w_next = w+1 then  	
		((W-w)*(1-defend_probab[a,b])/W)  
	else 0;					

var XIw {i in players, w in states} = 
	sum{a in actions} (if i=1 then XI1[w,a] else XI2[w,a])*sigma[i,w,a];


var PIw {w in states, w_ in states} = 
	sum{a1 in actions, a2 in actions} PI[w,w_,a1,a2]*sigma[1,w_,a1]*sigma[2,w_,a2];
	# Here too may be a limiation. DO recheck


minimize total_val:
	sum {i in players, w in states}(
		 XIw[i,w] + discount*sum{w_next in states}(PIw[w_next,w]*value[i,w_next]) - value[i,w] 
	);

subject to indiv_value {i in players, w in states,a1 in actions, a2 in actions}:
	(if i=1 then XI1[w,a1] else XI2[w,a2]) + discount*sum{w_next in states}(PI[w_next,w,a1,a2]*value[i,w_next]) >= value[i,w];
		
subject to probab_sum{i in players, w in states}:
	sum{a in actions} sigma[i,w,a]=1;

