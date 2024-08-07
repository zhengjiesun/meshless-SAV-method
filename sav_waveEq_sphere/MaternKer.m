function y = MaternKer(y,x,e,m,k)
global dim;


r = DistanceMatrix(y,x);
% k-th derivative of Matern kernel
r = r+eps;
switch k
    case 0
        y = 2^(1-(m-dim/2))/gamma(m-dim/2)*(e*r).^(m-dim/2).*besselk(m-dim/2,e*r);
    case 1
        y = -e^2*2^(1-(m-dim/2))/gamma(m-dim/2)*(e*r).^(m-1-dim/2).*besselk(m-1-dim/2,e*r);
end
