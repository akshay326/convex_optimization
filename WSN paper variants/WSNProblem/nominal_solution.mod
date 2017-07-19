param W;
param alpha;
param beta;
param inst_cost {w in 0..W} = alpha * (exp(beta*w)-1);

param Ch1;
param Ch2;
param Cl1;
param Cl2;

param discount;
set actions = {"low","high"};

param s1 {actions};
param s2 {actions};
param Pmin;
param Pmax;
param Pd = Uniform(Pmin,Pmax);
#var Pd <= Pmax, >= Pmin;


#IDS is represented by Player2
param chi1 {w in 0..W, a in actions, b in actions} = -inst_cost[w] + (if a="high" then Ch1 else Cl1);

param chi2 {w in 0..W, a in actions, b in actions} = inst_cost[w] + (if b="high" then Ch2 else Cl2);

param defend_probab {a in actions, b in actions} = ( if s1[a] =0 then 1 else min( 0.16*Pd*s2[b]/s1[a], 1)   );

param transition_probab {w in 0..W, w_next in 0..W, a in actions, b in actions} =  
					(if w_next = w-1 then (w*defend_probab[a,b]/W) 
					
					else (if w_next = w then (w*(1-defend_probab[a,b])/W  + (W-w)*defend_probab[a,b]/W)   
					
					else  (if w_next = w+1 then  ((W-w)*(1-defend_probab[a,b])/W)  else 0   )    ))  ;
					

var sigma {i in 1..2, w in 0..W, a in actions} >=0;

var C { i in 1..2, w in 0..W, a in actions, b in actions } = 
	( if i=1 then chi1[w,a,b]*sigma[2,w,b] else chi2[w,a,b]*sigma[1,w,a]  );

var transitn_probab_total {w in 0..W, w_next in 0..W} = sum { a in actions, b in actions } 
										(sigma[1,w,a]*sigma[2,w,b]* transition_probab[w,w_next,a,b] ) ;

var value {i in 1..2, w in 0..W} ;

minimize total_val : sum {i in 1..2, w in 0..W} 
						(  sum {a in actions, b in actions}
							 ( C[i,w,a,b]*(if i =1 then sigma[i,w,a] else sigma[i,w,b])  )  +
							 			discount* sum { w_next in 0..W} (transitn_probab_total[w, w_next] * value [i,w_next]) -value[i,w] );

subject to indiv_value {i in 1..2, w in 0..W,a in actions}:
value [i,w] <=  sum {b in actions} (if i =1 then C[i,w,a,b] else C[i,w,b,a] )   +
									discount* sum { w_next in 0..W} (
										(sum {b in actions} 
								( transition_probab[w, w_next,a,b]* (if i=1 then sigma[2,w,b] else sigma[1,w,b]) ) )* value [i,w_next] )   ;

		
subject to probab_sum { i in 1..2, w in 0..W } : sum {a in actions} sigma[i,w,a] = 1 ;

subject to nonnegative {i in 1..2, w in 0..W, a in actions} : sigma[i,w,a]>=0;


