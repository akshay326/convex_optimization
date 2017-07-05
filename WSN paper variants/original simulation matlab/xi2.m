function [x] = xi2( w,s2 )
%xi2 Returns cost for player 2(IDS)
%   Associated with current state 'w' and actions 's1','s2'
s2_low = 10;
s2_high = 50;

c2_low = 2;
c2_high = 10;

if s2 == s2_low
    x = -c(w) + c2_low;
else
    x = -c(w) + c2_high;
end
end

