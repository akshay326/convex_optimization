function [ value ] = pi_bar( w_bar,w,p )
%pi_bar Returns the state transition values
W = 30; % Size of state space

switch w_bar - w
    case -1
        value = w*p/W;
    case 0
        value = w*(1-p)/W + (1-(w/W))*p;
    case 1
        value = (1-(w/W))*(1-p);
    otherwise
        value = 0;
end

end

