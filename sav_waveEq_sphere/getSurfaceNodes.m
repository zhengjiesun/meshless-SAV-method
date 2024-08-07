function [X,Y,weight] = getSurfaceNodes(Nx,Ny)
global surfaceType
switch surfaceType
    case 'sphere'
        X = getMinEnergyNodes(Nx);
        % Y = getIcosNodes(Icosk(1),Icosk(2));Ny = size(Y,1);
        % %Y = X;Ny = size(Y,1);
        % %hX = 1/sqrt(Nx);
        % % For Y, we compute the quadrature weights, we use the kernel phi_m with m=2;
        % weight = sphere_quadweight(Y,2);

        switch Ny
            case 2562
                load('Icos2562.mat','weight','x');
            case 10242
                load('Icos10242.mat','weight','x');
            case 23042
                load('Icos23042.mat','weight','x');
            case 40962
                load('Icos40962.mat','weight','x');     
            case 2501
                load('sph2500.mat','weight_app','weight_exact','x');
            case 5001
                load('sph5000.mat','weight_app','weight_exact','x');
            case 7501
                load('sph7500.mat','weight_app','weight_exact','x');
            case 10001
                load('sph10000.mat','weight_app','weight_exact','x');
            case 15001
                load('sph15000.mat','weight_app','weight_exact','x');
            case 20001
                load('sph20000.mat','weight_app','weight_exact','x');
            case 30001
                load('sph30000.mat','weight_app','weight_exact','x');
            case 40001
                load('sph40000.mat','weight_app','weight_exact','x');  
            case 60001
                load('sph60000.mat','weight_app','weight_exact','x'); 
        end
        Y = x; 
        weight = weight_exact;
        
        
    case 'torus'
        switch Nx
            case 440
                load('torus440.mat','x'); X = x;
            case 724
                load('torus724.mat','x'); X = x;
            case 868
                load('torus868.mat','x'); X = x;
            case 1024
                load('torus1024.mat','x'); X = x;
            case 1380
                load('torus1380.mat','x'); X = x;
            case 1952
                load('torus1952.mat','x'); X = x;
            case 3112
                load('torus3112.mat','x'); X = x;
            case 3124
                load('torus3124.mat','x'); X = x;
            case 3968
                load('torus3968.mat','x'); X = x;
            case 5534
                load('torus5534.mat','x'); X = x;
            case 7520
                load('torus7520.mat','x'); X = x;
        end
        switch Ny
            case 1958
                load('torus1958.mat','weight_app','weight_exact','x');
            case 3124
                load('torus3124.mat','weight_app','weight_exact','x');
            case 5534
                load('torus5534.mat','weight_app','weight_exact','x');
            case 7520
                load('torus7520.mat','weight_app','weight_exact','x');
            case 9470
                load('torus9470.mat','weight_app','weight_exact','x');
             case 15560
                load('torus15560.mat','weight_app','weight_exact','x');
             case 21030
                load('torus21030.mat','weight_app','weight_exact','x');
            case 30704
                load('torus30704.mat','weight_app','weight_exact','x');
            case 39122
                load('torus39122.mat','weight_app','weight_exact','x');       
        end
     Y = x; weight = weight_exact;
    
     
     case 'ellipse'
        switch Nx
            case 1814
                load('ellipse1814.mat','x'); X = x;
        end
        switch Ny
            case 16230
                load('ellipse16230.mat','weight_app','weight_exact','x');
        end
     Y = x; weight = weight_exact;
     
     
    case 'cyclide'
        switch Nx
            case 1906
                load('cyclide1906.mat','x'); X = x;
           case 1634
                load('cyclide1634.mat','x'); X = x;
        end
        switch Ny
            case 7040
                load('cyclide7040.mat','weight_app','weight_exact','x');
            case 9818
                load('cyclide9818.mat','weight_app','weight_exact','x');
            case 11884
                load('cyclide11884.mat','weight_app','weight_exact','x');
        end
     Y = x; weight = weight_exact;
     
     
     case 'bretzel'
        switch Nx
            case 1648
                load('bretzel1648.mat','x'); X = x;
        end
        switch Ny
             case 8960
                load('bretzel8960.mat','weight_app','weight_exact','x');
            case 13014
                load('bretzel13014.mat','weight_app','weight_exact','x');
        end
     Y = x; weight = weight_exact;
     
    case 'rbc'
        X = getMinEnergyNodes(Nx);
        r0=3.39; c0=0.81/r0; c2=7.83/r0; c4=-4.39/r0; a=3.91/r0;
        rbc = @(la,th) [a*cos(la).*cos(th) a*sin(la).*cos(th) ...
            0.5*sin(th).*(c0+c2*(cos(th)).^2+c4*(cos(th)).^4)];
        % Use minimum energy nodes projected to the sphere
        [la,th]=cart2sph(X(:,1),X(:,2),X(:,3)); X=rbc(la,th);
        switch Ny
            case 2562
                load('Icos2562.mat','weight','x');
            case 10242
                load('Icos10242.mat','weight','x');
            case 23042
                load('Icos23042.mat','weight','x');
            case 40962
                load('Icos40962.mat','weight','x');     
            case 2501
                load('sph2500.mat','weight_app','weight_exact','x');
            case 5001
                load('sph5000.mat','weight_app','weight_exact','x');
            case 7501
                load('sph7500.mat','weight_app','weight_exact','x');
            case 10001
                load('sph10000.mat','weight_app','weight_exact','x');
            case 15001
                load('sph15000.mat','weight_app','weight_exact','x');
            case 20001
                load('sph20000.mat','weight_app','weight_exact','x');
            case 30001
                load('sph30000.mat','weight_app','weight_exact','x');
        end
        [la,th]=cart2sph(x(:,1),x(:,2),x(:,3)); Y=rbc(la,th);
        weight = weight_exact;
end
                    
    
end