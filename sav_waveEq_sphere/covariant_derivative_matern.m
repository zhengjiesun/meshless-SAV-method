function [A] = covariant_derivative_matern(y,x,ep,m)
% compute the covariant dderivative of matern kernel in Cartesian coordinates. 
% x, y could be one point or vectors.
% x center point,  y sample point
global surfaceType

N = size(x,1); M = size(y,1); sz = [M N];
DPx = zeros(sz); DPy = zeros(sz); DPz = zeros(sz); 


xi = x(:,1); yi = x(:,2); zi = x(:,3);
switch surfaceType
    case 'sphere'
        nxi = y(:,1);
        nyi = y(:,2);
        nzi = y(:,3);
        
    case 'torus' 
        nr = normalvec_torus(y,1,1/3);
        nxi = nr(:,1);nyi = nr(:,2);nzi = nr(:,3);

    case 'ellipse'
        a = 2; b = 1; c = 1.5;
        nxi = 2*y(:,1)/a^2;
        nyi = 2*y(:,2)/b^2;
        nzi = 2*y(:,3)/c^2;
        nr = sqrt(nxi.^2+nyi.^2+nzi.^2);
        nxi = nxi./nr; nyi = nyi./nr; nzi = nzi./nr;

    case 'cyclide'
        a = 2; b = 1.9; d = 1; c = sqrt(a^2-b^2);
        nxi = 4*y(:,1).*(sum(y.^2,2)-d^2+b^2)-8*a*(a*y(:,1)+c*d);
        nyi = 4*y(:,2).*(sum(y.^2,2)-d^2+b^2)-8*b^2*y(:,2);
        nzi = 4*y(:,3).*(sum(y.^2,2)-d^2+b^2);
        nr = sqrt(nxi.^2+nyi.^2+nzi.^2);
        nxi = nxi./nr; nyi = nyi./nr; nzi = nzi./nr;
        
    case 'bretzel'
        nxi = (4*y(:,1)-8*y(:,1).^3).*(y(:,1).^2.*(1-y(:,1).^2)-y(:,2).^2)-1/20*y(:,1);
        nyi = -4*y(:,2).*(y(:,1).^2.*(1-y(:,1).^2)-y(:,2).^2)-1/20*y(:,2);
        nzi = 19/20*y(:,3);
        nr = sqrt(nxi.^2+nyi.^2+nzi.^2);
        nxi = nxi./nr; nyi = nyi./nr; nzi = nzi./nr;
        
    case 'rbc'
        r0=3.39; c0=0.81/r0; c2=7.83/r0; c4=-4.39/r0; a=3.91/r0;
        rbc = @(la,th) [a*cos(la).*cos(th) a*sin(la).*cos(th) ...
            0.5*sin(th).*(c0+c2*(cos(th)).^2+c4*(cos(th)).^4)];
        % rbc = @(la,th) [(1+1/3*cos(th)).*cos(la) (1+1/3*cos(th)).*sin(la) (1+1/3*sin(th))];
        syms la th; % Compute the surface normal vectors symbolically
        nr = inline(cross(diff(rbc(la,th),la),diff(rbc(la,th),th)));
        % Use minimum energy nodes projected to the sphere
        [la,th]=cart2sph(y(:,1),y(:,2),y(:,3)); y=rbc(la,th);
        nr=nr(la,th); nr=nr./repmat(sqrt(sum(nr.^2,2)),[1 3]);
        nxi = nr(:,1);nyi = nr(:,2);nzi = nr(:,3);
        
end

D1M = MaternKer(y,x,ep,m,1);
 for k = 1:M
     DPx(k,:) = (1-nxi(k)^2)*(y(k,1)-xi) - nxi(k)*nyi(k)*(y(k,2)-yi) - nxi(k)*nzi(k)*(y(k,3)-zi);
     DPy(k,:) = -nxi(k)*nyi(k)*(y(k,1)-xi) + (1-nyi(k)^2)*(y(k,2)-yi) - nyi(k)*nzi(k)*(y(k,3)-zi);
     DPz(k,:) = -nxi(k)*nzi(k)*(y(k,1)-xi) - nyi(k)*nzi(k)*(y(k,2)-yi) + (1-nzi(k)^2)*(y(k,3)-zi);
 end

DPx =  D1M.*DPx; DPy = D1M.*DPy; DPz = D1M.*DPz;
A(:,:,1) = DPx; A(:,:,2) = DPy; A(:,:,3) = DPz;


end