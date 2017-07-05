function [x] = xi1( w,s1 )
%xi1 Returns cost for player 1(Intruder)
%   Associated with current state 'w' and actions 's1','s2'
s1_low = 0;
s1_high = 8;

c1_low = 0;
c1_high = 8;

if s1 == s1_low
    x = -c(w) + c1_low;
else
    x = -c(w) + c1_high;
end
end

