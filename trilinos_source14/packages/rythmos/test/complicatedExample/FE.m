%
% This is the exact output of Forward Euler applied to
% \dot{x} - \Lambda x = 0
% x(0) = x0
% t = 0 .. 1
% with fixed step-sizes equal to 1/n.
%
% Example:
% n = 4;
% lambda = [-2:3/(n-1):1]';
% x0 = 10*ones(n,1);
% FE(x0,lambda,n);
%
function xn = FE(x0,lambda,n)

xn = x0.*(1+lambda/n).^n;
