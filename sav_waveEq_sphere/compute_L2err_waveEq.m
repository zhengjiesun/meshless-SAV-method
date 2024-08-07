function err = compute_L2err_waveEq(Nx,sp,m,T, tau)
% This function aims to compute the L2err of the proposed method for
% solving wave equation on the sphere or the torus.

% Input  Nx: number of centers used in finite-dimensional approximation
%               sp: the shape parameter of the kernel
%               m: the degree of the kernel
%                T: the final time
%            tau: time step

%Output  err: L2error

        global surfaceType;

        switch surfaceType        
            case 'sphere'
                % number of evaulation point set
                Ny = 40001; 
                % exact solution
                exactu = @(t,lam,th) cos(t).*sin(th);
                % initial solution u0,v0
                testu = @(lam,th) sin(th); testv = @(lam,th) 0.*sin(th);
                %nonlinear term w.r.t solution
                Ffunc = @(x) -x.^2/2; d1Ffunc = @(x) -x; 
                % rhs nonlinear function
                nlf = @(x,y,z,t) 0*x.*y.*z;
                
            case 'torus'
                % number of evaulation point set
                Ny = 39122; 
                % exact solution
                exactu =  @(x,y,z,t)  cos(t)*x.*(x.^4-10*x.^2.*y.^2+5*y.^4)...
                    .*(x.^2+y.^2-60*z.^2);
                % initial solution u0,v0
                testu = @(x,y,z) x.*(x.^4-10*x.^2.*y.^2+5*y.^4)...
                    .*(x.^2+y.^2-60*z.^2);
                testv = @(x,y,z) 0*x.*y.*z;
                %nonlinear term w.r.t solution
                Ffunc = @(x) x.^2/2; d1Ffunc = @(x) x;
                % rhs nonlinear function
                nlf = @(x,y,z,t) -cos(t)*x.*(x.^4-10*x.^2.*y.^2+5*y.^4)...
                    .*(x.^2+y.^2-60*z.^2)+3*cos(t)./((sqrt(x.^2+y.^2)).^2).*...
                    (x.*(x.^4-10*x.^2.*y.^2+5*y.^4).*(10248*(sqrt(x.^2+y.^2)).^4-34335*(sqrt(x.^2+y.^2)).^3+...
                    41359*(sqrt(x.^2+y.^2)).^2-21320*(sqrt(x.^2+y.^2))+4000))+cos(t)*x.*(x.^4-10*x.^2.*y.^2+5*y.^4)...
                    .*(x.^2+y.^2-60*z.^2);        
        end
        
        
       %% Compute the mass matrix GramMat0, stiffness matrix GramMat1, Interpolation matrix IM
        [ctrs,epoints,GM0,GM1,IM,Func_lag,weight] = surface_matern_LLB(Nx,Ny,sp,m);

               
      %% Time-stepping under the framework of SAV by using Crank-Nicolson extrapolation
        Nt = ceil(T/tau);
        % define solution matrices and vectors
        Alp = zeros(Nx,Nt+1); Beta = zeros(Nx,Nt+1); p = zeros(Nt+1,1); 
        
        % compute the initial conditions and exact solution on evaulation points set
        switch surfaceType
            case 'sphere'
                [phix,thx] = cart2sph(ctrs(:,1),ctrs(:,2),ctrs(:,3));
                [Ephix,Ethx] = cart2sph(epoints(:,1),epoints(:,2),epoints(:,3));
                U0 = testu(phix,thx);V0 = testv(phix,thx); % Initial conditions
                uexact = exactu(T,Ephix,Ethx);
            case 'torus'
                U0 = testu(ctrs(:,1),ctrs(:,2),ctrs(:,3)); % Initial conditions
                V0 = testv(ctrs(:,1),ctrs(:,2),ctrs(:,3));
                uexact = exactu(epoints(:,1),epoints(:,2),epoints(:,3),T);
        end

       % transform the initial condition U0, V0 to its RBF expansion coefficients
        U0 = IM\U0; V0 = IM\V0;
        % new matrix in the algorithm
        Du = (GM0+tau^2/4*GM1); Dv = (GM0-tau^2/4*GM1);

        %  Compute the first time step
        C0 = 0;q0 =  sqrt(weight'*Ffunc(Func_lag*U0)+C0); p0 = q0;
        F0 = Func_lag'*(weight.*d1Ffunc(Func_lag*U0));
        gam0 = Dv*U0 + tau*GM0*V0 -tau^2/2*(p0/q0-U0'*F0/(4*q0^2))*F0...
            + tau^2/2*Func_lag'*(weight.*nlf(epoints(:,1),epoints(:,2),epoints(:,3),0));
        tempCoeff = Du\F0;
        tempf = gam0'*tempCoeff/(1+tau^2/(8*q0^2)*(F0'*tempCoeff));

        U1 = Du\(gam0-tau^2/(8*q0^2)*tempf*F0);
        V1 = 2*(U1-U0)/tau-V0;
        p1 = p0 + (U1-U0)'*F0/(2*q0);

        % save the initial solution and first time step solution in Alp, Beta, and p.
        Alp(:,1) = U0; Alp(:,2) = U1;
        Beta(:,1) = V0; Beta(:,2) = V1;
        p(1) = p0; p(2) = p1;
        
        
      %% Time-stepping using SAV approach
        for i = 2:Nt
            % extrapolation
            Alp_temp = 1/2*(3*Alp(:,i)-Alp(:,i-1)); 
            % compute the intermediate values
            q0 = sqrt(weight'*Ffunc(Func_lag*Alp_temp)+C0);
            F0 = Func_lag'*(weight.*d1Ffunc(Func_lag*Alp_temp));
            gam0 = Dv*Alp(:,i) + tau*GM0*Beta(:,i)- tau^2/2 * (p(i)/q0 -Alp(:,i)'*F0/(4*q0^2))*F0...
                + tau^2/2*Func_lag'*(weight.*nlf(epoints(:,1),epoints(:,2),epoints(:,3),(i-1/2)*tau));
            tempCoeff = Du\F0;
            tempf = (gam0'*tempCoeff)/(1+tau^2/(8*q0^2)*(F0'*tempCoeff));
            % update the solution at the next time level
            Alp(:,i+1) = Du\(gam0-tau^2/(8*q0^2)*tempf*F0);
            Beta(:,i+1) = 2*(Alp(:,i+1)-Alp(:,i))/tau-Beta(:,i);
            p(i+1) = p(i) + (Alp(:,i+1)-Alp(:,i))'*F0/(2*q0);
        end
        % compute the approximated solution at the evaluation point set
        uapp = Func_lag*Alp(:,end); 
        % compute the relative weighted L2 error
        err = sqrt(weight'*(uapp-uexact).^2)/sqrt(weight'*uexact.^2);
        
        
        %%  compute the discrete energy for each time level
%         for i = 1:Nt+1
%             Ener(i) = 1/2* (V(:,i)'*GramMat0*V(:,i) + U(:,i)'*GramMat1*U(:,i)) + p(i)^2;
%         end
%         plot(abs((Ener-Ener(1))/Ener(1)));
end