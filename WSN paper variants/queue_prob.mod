param S = 30;	# No of states
param I = 2;	# No of players
param J = 2;	# No of actions , assumed equal for all players

# 1 -> Router, 2 -> Service Provider
param lamda{i in 1..I} = 
	if  i = 1 then 1/10
	else 1/25;
	
param mu {i in 1..I} =
	if i =1 then Uniform(1/12,1/10)
	else Uniform(1/21,1/19);
	
param theta {1..J, 1..J};
let theta[1,1] := 90;
let theta[1,2] := 90;
let theta[2,1] := 110;
let theta[2,2] := 110;

param rho {1..J, 1..J};
let rho[1,1] := 70;
let rho[1,2] := 20;
let rho[2,1] := 30;
let rho[2,2] := 60;

param beta = 0.95;


/*
	Data Part Ends Here
	Model Starts
*/

var sigma {1..I, 0..S, 1..J};
var value {1..I, 0..S};


/*
	Creating PI and its weighted counterpart, PIw
	Note pg 11, Queue Paper, bottom right, eqn (18)
	I guess s=k=0 and s=k=S -> 1. The last conditions 
	written there are I guess wrong
*/

param PI {k in 0..S, s in 0..S, a1 in 1..J, a2 in 1..J} = 
	if s!=0 and k = s-1 then 
		mu[a1]/(lamda[a2] + mu[a1])  
	else if s != S and k = s+1 then
		lamda[a2]/(lamda[a2] + mu[a1])    
	else 0;


var PIw {w1 in 0..S, w2 in 0..S} = 
	sum {a in 1..J, b in 1..J} sigma[1,w1,a]*sigma[2,w1,b]*PI[w1,w2,a,b];


/* 
	Creating XI and its weighted counterpart, C[i,s,s(s)]
	
	Note:
	1) In WSN paper, Xi[i, s, s(s)] where s(s) is a vector of len J^I
		So, XI is I*S*J^I
	2) C[i,s,s(s)] following eqn (10) in WSN paper
		
*/

param h {s in 0..S} = 
	if s = 0 then 0
	else 1.2 * exp(0.2*s*log(1.9));

param XI {i in 1..I, s in 0..S, a1 in 1..J, a2 in 1..J } =	
	if i=1 then
		h[s] + theta[a1,a2]		# XI1 part
	else 
		rho[a1,a2];				# XI2 part

/*
	Weighted XI shud be matrix of J^I-1 rows n J cols		
	Ref #17, pg 9, top-right, eqn (8) states XI[a(-i), a(i)]
	The a(i) part can take values b/t 1..J [Columns]
	While a(-i) part takes J^(I-1) combinations [Rows]
	
	Here was the confusion!!
	In #17, this matrix is a function of 'i', 's'
	While in WSN paper, it depends on s(s) as well - which is wrong 

var C {i in 1..I, s in 0..S,
		a in 1..J,   	# If I=3, it had been a1 in 1..J, a2 in 1..J 
		col in 1..J } =	
		if i = 1 then
			XI[1,s,col,a]*sigma[2,s,col]
		else if i = 2 then
			XI[2,s,a,col]*sigma[1,s,col];
			
#		else if i = 3 then
#		XI[i,s,a1,a2,col]...
*/

var C {i in 1..I, s in 0..S, action in 1..J} =	
	if i = 1 then
		(XI[i,s,action,1] + XI[i,s,action,2])*sigma[2,s,action]
	else if i = 2 then
		(XI[i,s,1,action] + XI[i,s,2,action])*sigma[1,s,action];

/*
	minimize
	maximize
*/
minimize diff_to_ideal_value: 
		sum {i in 1..I, s in 0..S} (
			value[i,s]
			- (sum {l in 1..J} C[i,s,l]*sigma[i,s,l])
			- (beta*sum {w2 in 0..S} (value[i,w2]*PIw[w2, s]))
		);

subject to max_condition {i in 1..I, s in 0..S, a in 1..J}:
		value[i, s] >= C[i,s,a] 
			   + beta*sum {b in 1..J, w2 in 0..S} value[i,w2]*PI[s, w2,a,b]*(if i=1 then sigma[2,s,b] else sigma[1,s,b]);

subject to probab_sum {i in 1..I, s in 0..S}:
		sum {a in 1..J} (sigma[i,s,a]) = 1;
		
subject to positive_probab {i in 1..I, s in 0..S, a in 1..J} :
		sigma[i,s,a] >= 0;