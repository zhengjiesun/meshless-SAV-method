function [X,Y,A,B,IM,func_lag,weight] = surface_matern_LLB(Nx,Ny,ep,m)
%Input      Nx : number of centers used for constructing RBF approximation
%                 Ny:  number of quadrature points used for spherical integral
%                 ep:   shape parameter of RBF kernel,  m: order of kernel

%Output    X: centers    Y: quadrature points
%                  A: mass matrix G_0;
%                  B: stiffness matrix G_1;
%                  IM: interpolation matrix
%                  func_lag: evaulation matrix 
%                  weight: quadrature weight

%% Node set, X for interpolation, Y for quadrature, weight for quadrature rule
[X,Y,weight] = getSurfaceNodes(Nx,Ny);

%% Compute coefficients of local Lagrange function
M = 3;  n = M*ceil((log10(Nx))^2);
Md1 = createns(X,'Distance','cosine');
idx = knnsearch(Md1,X,'K',n,'Distance','cosine');
% compute local weights
W = zeros(Nx,n);
for k = 1:Nx
    W(k,:) = compute_weights_matern(X(idx(k,:),:),n,ep,m);
end

KX = MaternKer(X,X,ep,m,0); 
IM = zeros(Nx,Nx);
% compute the Lagrange function value at Y for each X
KY = MaternKer(Y,X,ep,m,0); 
func_lag = zeros(Ny,Nx);
for k = 1:Nx
    IM(:,k) = KX(:,idx(k,:))*W(k,:)';
    func_lag(:,k) = KY(:,idx(k,:))*W(k,:)'; % Lagrange function value
end
%-------------------------------------------------------------------------

%% compute the covariant derivatives
[Dphi] = covariant_derivative_matern(Y,X,ep,m); 
Dx = Dphi(:,:,1); Dy = Dphi(:,:,2); Dz = Dphi(:,:,3);

for k = 1:Nx
    DDx(:,k) = Dx(:,idx(k,:))*W(k,:)';
    DDy(:,k) = Dy(:,idx(k,:))*W(k,:)';
    DDz(:,k) = Dz(:,idx(k,:))*W(k,:)';
end
Dx = DDx; Dy = DDy; Dz = DDz;
% Define the stiffness matrix
A = zeros(Nx,Nx); B = zeros(Nx,Nx);
for j = 1:Nx
    for k = j:Nx
        gradient_term = Dx(:,j).*Dx(:,k)+Dy(:,j).*Dy(:,k)+Dz(:,j).*Dz(:,k);
        func_term = func_lag(:,j).*func_lag(:,k);
        A(j,k) = func_term'*weight;
        B(j,k) = gradient_term'*weight;
    end
end
A = (A + A')-diag(diag(A)); B = (B + B')-diag(diag(B));
