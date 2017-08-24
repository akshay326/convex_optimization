# option minos_options 'meminc=15' 'iterations_limit=100';
# option loqo_solver 'iterlim=100';
# option solver LOQO;
# option solver MINOS;

# reset;option solver LOQO;model 25_8_17.mod;data 25_8_17.dat;solve;
/*
	:( -> Increases dual objective
	:) -> Decreases dual objective
	:| -> No diff
	
	primal objective 13.66772521
  	dual objective -52,212.02834
	
	Possible Edits: 
		
		
		
		
	Edits Evaluated:
		
		Commenting definition of Y -> :(
		Turn instance 1 to 2 -> :(
		Turn e4 to objective - :)
		Commenting e3 -> :|
		Turn value[i,w] to value[w], *-1 for router -> :(
		Turn beta*Y*x to beta*z in In3 -> :(
		
*/

param S;
param i_n = 1;
param beta;
 
set players;
set states = {0..30};
set instances = {1..5};
set actions;
set lim_type;

param lamda{actions};	# SrvcPro -> 1 Router  -> 2
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
		
		

# Refer pg 12, Queue Paper. Theta is indexed as {mu,lamda}
param theta {actions, actions};
param rho {actions, actions};

/*
	Fig 2, pg 13: for{w in states}{print(sum{a in actions} lamda[a]*sigma['Router',w,a]);}
*/
var sigma {players, states, actions} >= 0;
var value {players,states};


param PILower {k in states,s in states, a1 in actions, a2 in actions} =  
	if s=0 and k=1 then
		1
	else if s=S and k=S-1 then
		1
	else if s!=0 and k = s-1 then 
		mu[a1,i_n,'lower']/(lamda[a2] + mu[a1,i_n,'lower'])  
	else if s!=S and k = s+1 then
		lamda[a2]/(lamda[a2] + mu[a1,i_n,'upper']);
			
param PIUpper {k in states,s in states, a1 in actions, a2 in actions} =
	if s=0 and k=1 then
		1
	else if s=S and k=S-1 then
		1
	else if s!=0 and k = s-1 then 
		mu[a1,i_n,'upper']/(lamda[a2] + mu[a1,i_n,'upper'])  
	else if s!=S and k = s+1 then
		lamda[a2]/(lamda[a2] + mu[a1,i_n,'lower']);

param h {s in states} = 
	if s = 0 then 0
	else 1.2 * exp(0.2*s*log(1.9));

# Profit to SrvcPro -> 1
param XI{s in states, a1 in actions, a2 in actions } = 
	h[s] + theta[a1,a2] - rho[a1,a2];


/*
	For debugging
 	display t.lb,t,t.ub,t.rc;
*/
var t{i in players,k in states, w in states, a1 in actions, a2 in actions} >= PILower[k,w,a1,a2], <= PIUpper[k,w,a1,a2];
var T{i in players, k in states, s in states}=
	sum{a1 in actions, a2 in actions} sigma['SrvcPro',s,a1]*sigma['Router',s,a2]*t[i,k,s,a1,a2];

var E {i in players, w in states, otherPlayer in actions, currPlayer in actions} = 
	if i = 'SrvcPro' then
		XI[w,currPlayer,otherPlayer]*sigma['Router',w,otherPlayer]
	else
		-XI[w,otherPlayer,currPlayer]*sigma['SrvcPro',w,otherPlayer];





param b{lim in lim_type, k in states,s in states, a1 in actions, a2 in actions} =
	(if lim = 'upper' then -PIUpper[k,s,a1,a2] else PILower[k,s,a1,a2]);
	
param A{w in states, lim in lim_type,
		w1 in states, a11 in actions, a12 in actions,
		w2 in states, a21 in actions, a22 in actions}=
	if a11 = a21 and a12 = a22 and w1=w2 then
		if lim = 'lower' then 
			1
		else
			-1
	else 0;
	
/*
	Just a small check
	print(sum {lim in lim_type,w in states,w1 in states, a1 in actions, a2 in actions} b[lim,w,w1,a1,a2]);
	print(sum {lim in lim_type,w in states,w1 in states, a1 in actions, a2 in actions}(sum{w2 in states, s1 in actions, s2 in actions} A[w,lim,w1,a1,a2,w2,s1,s2]*PI[lim,w,w2,s1,s2]));
	They must be same
*/
			

param Q{w in states, a1 in actions, a2 in actions, 
	w_next in states, s1 in actions, s2 in actions} = 
	if a1 = s1 and a2 = s2 then 
		1
	else 
		0;

var z{i in players, s in states, k in states, a1 in actions, a2 in actions}=
	value[i,k]*sigma['SrvcPro',s,a1]*sigma['Router',s,a2];
	
var Y{i in players, s in states, k in states, a1 in actions, a2 in actions, l in actions}=
	z[i,s,k,a1,a2]*sigma[i,s,l]/(1-(2*sigma[i,s,'low']*sigma[i,s,'high']));


var m{players, states, lim_type, states, actions, actions} <= 0;
var n{players, states, actions, actions};








minimize nothing:
	sum{i in players,s in states}
	(value[i,s] - (sum{other in actions, curr in actions}(E[i,s,other,curr]*sigma[i,s,curr]) + beta*sum{aux in states}(T[i,s,aux]*value[i,aux])));

subject to in1 {i in players,curr in actions, w in states}:
	  (sum{other in actions} E[i,w,other,curr]) 
	+ (beta*sum {w_next in states, s1 in actions, s2 in actions} Y[i,w,w_next,s1,s2,curr]*t[i,w,w_next,s1,s2])    
	>= value[i,w];


subject to in2 {i in players,w in states}:
	value[i,w] >= 
	  (sum {other in actions, curr in actions} E[i,w,other,curr]*sigma[i,w,curr])
	+ (sum{lim in lim_type, w_next in states, s1 in actions, s2 in actions} b[lim,w,w_next,s1,s2]*m[i,w,lim,w_next,s1,s2])
	+ (sum{s1 in actions, s2 in actions} n[i,w,s1,s2]);

	

subject to in3 {i in players, w in states, w1 in states, a1 in actions, a2 in actions}:
	(sum{lim in lim_type, w2 in states, s1 in actions, s2 in actions} m[i,w,lim,w2,s1,s2]*A[w,lim,w2,s1,s2,w1,a1,a2])
   +(sum{s1 in actions, s2 in actions} n[i,w,s1,s2]*Q[w,s1,s2,w1,a1,a2])
	>= beta*(Y[i,w,w1,a1,a2,'low']*sigma[i,w,'low'] + Y[i,w,w1,a1,a2,'high']*sigma[i,w,'high']); 





subject to e1 {i in players, s in states}:
	sigma[i,s,'low'] + sigma[i,s,'high'] = 1; 
	
subject to e2 {i in players, w in states, a1 in actions, a2 in actions} :
	(sum {w_next in states, s1 in actions, s2 in actions} Q[w,a1,a2,w_next,s1,s2]*t[i,w_next,w,s1,s2]) = 1; 
	
#subject to e3 {i in players, s in states, k in states, a1 in actions, a2 in actions}:
#	(sum{l in actions}Y[i,s,k,a1,a2,l]*sigma[i,s,l])= z[i,s,k,a1,a2];
	
subject to e4 {i in players,s in states}:
	value[i,s] >= sum{other in actions, curr in actions}(E[i,s,other,curr]*sigma[i,s,curr]) + beta*sum{aux in states}(T[i,s,aux]*value[i,aux]);
	