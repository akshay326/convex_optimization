function [ x ] = c(w)
%c Function used as penalty in Original Simulation
%  Uses Exponential penalty
%  'w' is the state
alpha = 10;
beta = 0.05;
x = alpha*(exp(beta*w) - 1);
end

