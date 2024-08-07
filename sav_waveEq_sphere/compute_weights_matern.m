function [w,rhs,CA] = compute_weights_matern(x,n,ep,m)

% Input x : local points 
%        n : number of local points
% use the matern kernel;
K = MaternKer(x,x,ep,m,0);  
rhs = [1;zeros(n-1,1)];
w = K\rhs;
CA = cond(K);