# reset;model 18_8_17.mod;data 18_8_17.dat;solve;

param S;
param beta;
 
set players;
set states = {0..30};
set actions;
set type;

param mu{actions};		# SrvcPro -> 1
param lamda{actions};	# Router  -> 2

# Refer pg 12, Queue Paper. Theta is indexed as {mu,lamda}
param theta {actions, actions};
param rho {actions, actions};

/*
	Fig 2, pg 13: for{w in states}{print(sum{a in actions} lamda[a]*sigma['Router',w,a]);}
	Fig 3a		: for{w in states}{print(sum{a in actions} mu[a]*sum{w_ in states}(sigma['SrvcPro',w_,a]*T['SrvcPro',w_,w]));}
	Fig 3b		: for{w in states}{print(sum{a in actions} mu[a]*sum{w_ in states}(sigma['SrvcPro',w_,a]*T['Router',w_,w]));}
*/
var sigma {players, states, actions} >= 0;
var value {players, states};


param PI {k in states,s in states, a1 in actions, a2 in actions} = 
	if s=0 and k=1 then
		1
	else if s=S and k=S-1 then
		1
	else if s!=0 and k = s-1 then 
		mu[a1]/(lamda[a2] + mu[a1])  
	else if s!=S and k = s+1 then
		lamda[a2]/(lamda[a2] + mu[a1]);

param h {s in states} = 
	if s = 0 then 0
	else 1.2 * exp(0.2*s*log(1.9));

# Profit to SrvcPro -> 1
param XI{s in states, a1 in actions, a2 in actions } = 
	h[s] + theta[a1,a2] - rho[a1,a2];

var t{i in players,k in states, w in states, a1 in actions, a2 in actions}=
	if i='SrvcPro' then 
		if (XI[w,'low',a2] + beta*sum{w_ in states} value['SrvcPro',w_]*PI[w_,w,'low',a2]) > (XI[w,'high',a2] + beta*sum{w_ in states} value['SrvcPro',w_]*PI[w_,w,'high',a2]) then
			PI[k,w,'low',a2]
		else
			PI[k,w,'high',a2]
	else
		if (-XI[w,a1,'low'] + beta*sum{w_ in states} value['Router',w_]*PI[w_,w,a1,'low']) > (-XI[w,a1,'high'] + beta*sum{w_ in states} value['Router',w_]*PI[w_,w,a1,'high']) then
			PI[k,w,a1,'low']
		else
			PI[k,w,a1,'high'];

var T{i in players, k in states, s in states}=
	sum{a1 in actions, a2 in actions} sigma['SrvcPro',s,a1]*sigma['Router',s,a2]*t[i,k,s,a1,a2];

var E {i in players, w in states, a1 in actions, a2 in actions} = 
	if i = 'SrvcPro' then
		XI[w,a1,a2]*sigma['Router',w,a2]
	else
		-XI[w,a1,a2]*sigma['SrvcPro',w,a1];


param b{_t in type,k in states,s in states, a1 in actions, a2 in actions}= 
	if _t = 'lower' then  
		if s=0 and k=1 then
			-1
		else if s=S and k=S-1 then
			-1
		else if s!=0 and k = s-1 then 
			-mu[a1]/(lamda[a2] + mu[a1])  
		else if s!=S and k = s+1 then
			-lamda[a2]/(lamda[a2] + mu[a1])
	else
		if s=0 and k=1 then
			1
		else if s=S and k=S-1 then
			1
		else if s!=0 and k = s-1 then 
			mu[a1]/(lamda[a2] + mu[a1])  
		else if s!=S and k = s+1 then
			lamda[a2]/(lamda[a2] + mu[a1]);


param A{w in states, _t in type, 
		w1 in states, a11 in actions, a12 in actions,
		w2 in states, a21 in actions, a22 in actions}=
	if a11 = a21 and a12 = a22 and w1=w2 then 
		if b[_t,w,w1,a11,a12] < 0  then  
			-1
		else if b[_t,w,w1,a11,a12] > 0  then
			1
		else
			0
	else 0;
			

#	Q(w) is L^N x M*L^N. So Q shud be M*L^N x M*L^N such that
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


var m{players, states, actions, states, actions, actions} <= 0;
var n{players, states, actions, actions};


minimize obj: 
	sum{i in players,w in states} value[i,w];

subject to in1 {w in states, a1 in actions}:
	value['SrvcPro',w] >= sum{a2 in actions} sigma['Router',w,a2]*(XI[w,a1,a2] + beta*sum{w_ in states} value['SrvcPro',w_]*PI[w_,w,a1,a2]);

subject to in2 {w in states, a2 in actions}:
	value['Router',w] >= sum{a1 in actions} sigma['SrvcPro',w,a1]*(-XI[w,a1,a2] + beta*sum{w_ in states} value['Router',w_]*PI[w_,w,a1,a2]);

/*
var shit;
for{i in players, w in states, _t in type, w1 in states, a1 in actions, a2 in actions}{
	let shit := (sum{w2 in states, s1 in actions, s2 in actions} A[w,_t,w1,a1,a2,w2,s1,s2]*t[i,w,w2,s1,s2])-b[_t,w,w1,a1,a2]);
	if shit!=0 then 
		display shit;
}
*/
subject to in6 {i in players, w in states, _t in type, w1 in states, a1 in actions, a2 in actions} :
	(sum{w2 in states, s1 in actions, s2 in actions} A[w,_t,w1,a1,a2,w2,s1,s2]*t[i,w,w2,s1,s2]) >=	b[_t,w,w1,a1,a2];


subject to e1 {i in players, s in states}:
	sigma[i,s,'low'] + sigma[i,s,'high'] = 1; 
	
subject to e2 {i in players, w in states, a1 in actions, a2 in actions} :
	(sum {w_next in states, s1 in actions, s2 in actions} Q[w,a1,a2,w_next,s1,s2]*t[i,w_next,w,s1,s2]) = 1; 
	
subject to e3 {i in players, s in states, k in states, a1 in actions, a2 in actions}:
	(sum{l in actions}Y[i,s,k,a1,a2,l]*sigma[i,s,l])= z[i,s,k,a1,a2];
	