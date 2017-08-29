param S;
param i_n = 1;
param beta;
 
set players;
set states = {0..30};
set instances = {1..5};
set actions;
set lim_type;

param lamda{actions};
param theta {actions, actions};
param rho {actions, actions};

param mu{a in actions, i in instances, lim in lim_type}=
	if a = 'low' then 
		if lim = 'lower' then
			1/(19 + 2*i)
		else
			1/(21 - 2*i)
	else
		if lim = 'lower' then
			1/(10 + 2*i)
		else
			1/(12 - 2*i);
	


param PILower {w in states, a1 in actions, a2 in actions,w_next in states} =  
	if w=0 and w_next=1 then
		1
	else if w=S and w_next=S-1 then
		1
	else if w!=0 and w_next = w-1 then 
		(1-(lamda[a2]/(lamda[a2] + mu[a1,i_n,'lower'])))  
	else if w!=S and w_next = w+1 then
		lamda[a2]/(lamda[a2] + mu[a1,i_n,'upper']);
			
param PIUpper {w in states, a1 in actions, a2 in actions,w_next in states} =
	if w=0 and w_next=1 then
		1
	else if w=S and w_next=S-1 then
		1
	else if w!=0 and w_next = w-1 then 
		(1-(lamda[a2]/(lamda[a2] + mu[a1,i_n,'upper'])))  
	else if w!=S and w_next = w+1 then
		lamda[a2]/(lamda[a2] + mu[a1,i_n,'lower']);

param h {w in states} = 
	if w = 0 then 0
	else 1.2 * exp(0.2*w*log(1.9));

# Profit to SrvcPro
param XI{w in states, a1 in actions, a2 in actions } = 
	h[w] + theta[a1,a2] - rho[a1,a2];




var sigma {players, states, actions};
var value {players,states};
var t{i in players, w in states, a1 in actions, a2 in actions,w_next in states} >= PILower[w,a1,a2,w_next], <= PIUpper[w,a1,a2,w_next];
var T{i in players, w in states, w_next in states}=
	sum{a1 in actions, a2 in actions} sigma['SrvcPro',w,a1]*sigma['Router',w,a2]*t[i,w,a1,a2,w_next];

var E {i in players, w in states, otherPlayer in actions, currPlayer in actions} = 
	if i = 'SrvcPro' then
		XI[w,currPlayer,otherPlayer]*sigma['Router',w,otherPlayer]
	else
		-XI[w,otherPlayer,currPlayer]*sigma['SrvcPro',w,otherPlayer];



param A{w in states, lim in lim_type,
		a11 in actions, a12 in actions, w1 in states,
		a21 in actions, a22 in actions, w2 in states}=
	if a11 = a21 and a12 = a22 and w1=w2 then
		if lim = 'lower' then 
			1
		else
			-1
	else 0;

param b{w in states, lim in lim_type, 
		a1 in actions, a2 in actions, w_next in states} =
	(if lim = 'upper' then -PIUpper[w,a1,a2,w_next] else PILower[w,a1,a2,w_next]);
	
			

param Q{w in states, 
		a1 in actions, a2 in actions, 
		s1 in actions, s2 in actions, w_next in states} = 
	if a1 = s1 and a2 = s2 then 
		1
	else 
		0;

var z{i in players, w in states, 
	  a1 in actions, a2 in actions, w_next in states}=
	value[i,w_next]*sigma['SrvcPro',w,a1]*sigma['Router',w,a2];
	
var Y{i in players, w in states,
	  a1 in actions, a2 in actions, w_next in states,
	  l in actions}=
	z[i,w,a1,a2,w_next]*sigma[i,w,l]/(1-(2*sigma[i,w,'low']*sigma[i,w,'high']));


var m{players, states, lim_type, actions, actions, states} <= 0;
var n{players, states, actions, actions};








minimize nothing:0;

subject to in1 {i in players,curr in actions, w in states}:
	  (sum{other in actions} E[i,w,other,curr]) 
	+ (beta*sum {w_next in states, s1 in actions, s2 in actions} Y[i,w,s1,s2,w_next,curr]*t[i,w,s1,s2,w_next])    
	>= value[i,w];


subject to in2 {i in players,w in states}:
	value[i,w] >= 
	  (sum {other in actions, curr in actions} E[i,w,other,curr]*sigma[i,w,curr])
	+ (sum{lim in lim_type, w_next in states, s1 in actions, s2 in actions} b[w,lim,s1,s2,w_next]*m[i,w,lim,s1,s2,w_next])
	+ (sum{s1 in actions, s2 in actions} n[i,w,s1,s2]);

	
subject to in3 {i in players, w in states, a1 in actions, a2 in actions, w_next in states}:
	(sum{lim in lim_type, s1 in actions, s2 in actions, w2 in states} A[w,lim,a1,a2,w_next,s1,s2,w2]*m[i,w,lim,s1,s2,w2])
   +(sum{s1 in actions, s2 in actions} Q[w,a1,a2,s1,s2,w_next]*n[i,w,s1,s2])
	>= beta*z[i,w,a1,a2,w_next];

subject to in4{i in players, w in states, a in actions}:
	1 >= sigma[i,w,a] >= 0;



subject to e1 {i in players, w in states}:
	sigma[i,w,'low'] + sigma[i,w,'high'] = 1; 

	
subject to e2 {i in players, w in states, a1 in actions, a2 in actions} :
	(sum {s1 in actions, s2 in actions, w_next in states} Q[w,a1,a2,s1,s2,w_next]*t[i,w,s1,s2,w_next]) = 1; 

#subject to e3 {i in players, w in states, w_next in states, a1 in actions, a2 in actions}:
#		(sum{l in actions}Y[i,w,w_next,a1,a2,l]*sigma[i,w,l])= z[i,w,w_next,a1,a2];

subject to e4{i in players,w in states}:
	 (if i='Router' then 10 else 0) >= value[i,w] - (sum{other in actions, curr in actions}(E[i,w,other,curr]*sigma[i,w,curr]) + beta*sum{aux in states}(T[i,w,aux]*value[i,aux])) >= (if i='Router' then 0 else -10);	