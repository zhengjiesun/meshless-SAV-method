
addpath nodes;
addpath sav_waveEq_sphere;

clear;
global  dim; dim = 3;     % dimension of the problem
global surfaceType;  
surfaceType = 'torus';  % surfaceType = sphere or torus


% Define problem parameters
mvec = [3, 4, 5]; % Order of Matern kernel
switch surfaceType
    case 'sphere'
        Nxvec = [21^2, 26^2, 32^2, 45^2, 56^2];
        spvec = [4, 6, 10];
    case 'torus'
        Nxvec = [440, 724, 1024, 1952, 3112];
        spvec = [3, 8, 8];
end

lm = length(mvec);
lx = length(Nxvec);
L2err = zeros(lx, lm);

% Calculate L2 error for each grid size and Matern kernel order
for j = 1:lm
    m = mvec(j) - 1;
    sp = spvec(j);
    for k = 1:lx
        Nx = Nxvec(k);
        L2err(k, j) = compute_L2err_waveEq(Nx, sp, m);
        fprintf('L2-err of wave equation Nx(%d) is: %e\n', Nx, L2err(k, j));
    end
end


