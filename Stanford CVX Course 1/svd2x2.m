%% Matlab function for computing the VSD of a 2x2 matrix
%% Written by Randy Ellis
%% More information at : http://www.lucidarme.me


function [U,SIG,V] = svd2x2(A)
% [U,SIG,V] = svd2x2(A) finds the SVD of 2x2 matrix A
% where U and V are orthogonal, SIG is diagonal,
% and A=U*SIG*V’
% Find U such that U*A*A’*U’=diag
Su = A*A';
phi = 0.5*atan2(Su(1,2)+Su(2,1), Su(1,1)-Su(2,2));
Cphi = cos(phi);
Sphi = sin(phi);
U = [Cphi -Sphi ; Sphi Cphi];


% Find W such that W’*A’*A*W=diag
Sw = A'*A;
theta = 0.5*atan2(Sw(1,2)+Sw(2,1), Sw(1,1)-Sw(2,2));
Ctheta = cos(theta);
Stheta = sin(theta);
W = [Ctheta, -Stheta ; Stheta, Ctheta];


% Find the singular values from U
SUsum= Su(1,1)+Su(2,2);
SUdif= sqrt((Su(1,1)-Su(2,2))^2 + 4*Su(1,2)*Su(2,1));
svals= [sqrt((SUsum+SUdif)/2) sqrt((SUsum-SUdif)/2)];
SIG =diag(svals);


% Find the correction matrix for the right side
S = U'*A*W
C = diag([sign(S(1,1)) sign(S(2,2))]);
V = W*C;