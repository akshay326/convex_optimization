function [ probab ] = p( s1,s2 )
%p Returns probability associated with s1, s2 actions
%   Currently taking uniform distribution
p_low = 0.8;
p_high = 0.8;
PD = unifrnd(p_low,p_high);

probab = min((PD*s2)/(50*s1));
end

