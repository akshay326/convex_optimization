# option minos_options 'meminc=15' 'iterations_limit=100';
# option loqo_solver 'iterlim=100';
# option solver LOQO/MINOS;

# reset;model 22_8_17.mod;data 22_8_17.dat;solve;

# Today i'll try creating theorem 5

param S;
param beta;
param EPS = 0.0001;
 
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
var value {players, states};


param PILower {k in states,s in states, a1 in actions, a2 in actions} =  
	if s=0 and k=1 then
		1
	else if s=S and k=S-1 then
		1
	else if s!=0 and k = s-1 then 
		mu[a1,1,'lower']/(lamda[a2] + mu[a1,1,'upper'])  
	else if s!=S and k = s+1 then
		lamda[a2]/(lamda[a2] + mu[a1,1,'upper']);
			
param PIUpper {k in states,s in states, a1 in actions, a2 in actions} =
	if s=0 and k=1 then
		1
	else if s=S and k=S-1 then
		1
	else if s!=0 and k = s-1 then 
		mu[a1,1,'upper']/(lamda[a2] + mu[a1,1,'lower'])  
	else if s!=S and k = s+1 then
		lamda[a2]/(lamda[a2] + mu[a1,1,'lower']);

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

var E {i in players, w in states, a1 in actions, a2 in actions} = 
	if i = 'SrvcPro' then
		XI[w,a1,a2]*sigma['Router',w,a2]
	else
		-XI[w,a1,a2]*sigma['SrvcPro',w,a1];





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
	
# If i dont define it n put it in constrains, it takes a lot of time
# Y[..L] = z[..]*u[L]/(u[L]^2 + u[H]^2)
var Y{i in players, s in states, k in states, a1 in actions, a2 in actions, l in actions}=
		z[i,s,k,a1,a2]*sigma[i,s,l]/(1-(2*sigma[i,s,'low']*sigma[i,s,'high']));


var m{players, states, lim_type, states, actions, actions} <= 0;
var n{players, states, actions, actions};








minimize nothing:0;

subject to in1_SrvcPro {a1 in actions, w in states}:
	  (sum{a2 in actions} E['SrvcPro',w,a1,a2]) 
	+ (beta*sum {w_next in states, s1 in actions, s2 in actions} Y['SrvcPro',w,w_next,s1,s2,a1]*t['SrvcPro',w,w_next,s1,s2])    
	>= value['SrvcPro', w];


subject to in1_Router {a2 in actions, w in states}:
	  (sum{a1 in actions} E['Router',w,a1,a2]) 
	+ (beta*sum {w_next in states, s1 in actions, s2 in actions} Y['Router',w,w_next,s1,s2,a2]*t['Router',w,w_next,s1,s2])    
	>= value['Router', w];
	

subject to in2 {i in players, w in states}:
	value[i, w] >= 
	  (sum {s1 in actions, s2 in actions} E[i,w,s1,s2]*(if i ='SrvcPro' then sigma[i,w,s1] else sigma[i,w,s2]) )
	+ (sum{lim in lim_type, w_next in states, s1 in actions, s2 in actions} b[lim,w,w_next,s1,s2]*m[i,w,lim,w_next,s1,s2])
	+ (sum{s1 in actions, s2 in actions} n[i,w,s1,s2]);

subject to in3 {i in players, w in states, w1 in states, a1 in actions, a2 in actions}:
	(sum{lim in lim_type, w2 in states, s1 in actions, s2 in actions} m[i,w,lim,w2,s1,s2]*A[w,lim,w2,s1,s2,w1,a1,a2])
   +(sum{s1 in actions, s2 in actions} n[i,w,s1,s2]*Q[w,s1,s2,w1,a1,a2])
	>= beta*z[i,w,w1,a1,a2]; 


subject to in6 {i in players,lim in lim_type,w in states,w1 in states, a1 in actions, a2 in actions} : 
	(sum{w2 in states, s1 in actions, s2 in actions} A[w,lim,w1,a1,a2,w2,s1,s2]*t[i,w,w2,s1,s2]) >=	b[lim,w,w1,a1,a2];





subject to e1 {i in players, s in states}:
	sigma[i,s,'low'] + sigma[i,s,'high'] = 1; 
	
subject to e2 {i in players, w in states, a1 in actions, a2 in actions} :
	(sum {w_next in states, s1 in actions, s2 in actions} Q[w,a1,a2,w_next,s1,s2]*t[i,w_next,w,s1,s2]) = 1; 
	
subject to e3 {i in players, s in states, k in states, a1 in actions, a2 in actions}:
	(sum{l in actions}Y[i,s,k,a1,a2,l]*sigma[i,s,l])= z[i,s,k,a1,a2];
	
subject to e4 {i in players, s in states}:
	value[i,s] = sum{row in actions, col in actions}(E[i,s,row,col]*sigma[i,s,row]) + beta*sum{aux in states}(T[i,s,aux]*value[i,aux]);