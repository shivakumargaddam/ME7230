classdef Solver
    % HELP - Solver (Obtain Solution)
    %   -- Methods:
    %               stagPFM [the staggered algorithm]
    %               elasticPFM [solves the elastic problem]
    %               damagePFM [solves the phase field problem]
    %               Psi_pm [calculates strain and strain energy]
    %               getShape2d [shape functions for 2d element: triangular,quad4,quad8]
    %               gauss_pw2D [generates Gauss Points for 2D analysis]
    %               lgwt [Legendre Gauss Weights in 1D]

    properties
    end

    methods (Static)
        %
        function stagPFM(duincMAX,duincMIN,duinc)
            % Extraction of Variables -------------------------------------
            simdata         = matfile('tempresults.mat','Writable',true);
            eType           = simdata.eType;
            nodes           = simdata.nodes;
            elements        = simdata.elements;
            uebcnodes       = simdata.uebcnodes;
            uebcvalues      = simdata.uebcvalues;
            phiebcnodes     = simdata.phiebcnodes;
            phiebcvalues    = simdata.phiebcvalues;
            ntsteps         = simdata.ntsteps;
            Charlen         = simdata.Charlen;
            EnergyGc        = simdata.EnergyGc;
            RatPoisson      = simdata.RatPoisson;
            YoungsMod       = simdata.YoungsMod;
            tolerance       = simdata.tolerance;
            reacnodes       = simdata.reacnodes;
            cracknodes      = simdata.cracknodes;
            crackvalues     = simdata.crackvalues;

            % Declaration of Field Variables ------------------------------
            noN = size(nodes,1);                                            % No of nodes
            noE = size(elements,1);                                         % number of elements
            uvdispevol(1:2*noN,1:ntsteps+1) = NaN;                          % Defining Global displacement vector
            phievol(1:noN,1:ntsteps+1) = NaN;                               % Defining Global phase field variable
            strain(1:3,1:noE,1:ntsteps+1) = NaN;                            % Defining strain for each element
            ReactionF(1:length(reacnodes),1:ntsteps+1) = NaN;               % Defining reaction force at fixed nodes
            Hplus(1:noE,1) = 0;                                             % History Variable
            Disp = zeros(ntsteps+1,1);                                      % Displacement (X-axis)

            % Initialization with an unstretched and undamaged state ------
            istep = 0;                                                      % at time t=0
            uvdispevol(:,istep+1) = 0;                                      % initial displacement
            phievol(:,istep+1) = 0;                                         % initial phase field
            phievol(cracknodes,istep+1) = crackvalues;                      % embedding initial cracks
            strain(:,:,istep+1) = 0;                                        % initial strain field
            ReactionF(:,istep+1) = 0;                                       % Defining reaction force at fixed nodes
            Disp(istep+1) = 0;                                              % initial displacement (X-axis in displacement vs. force plot)
            damage_switch = 0;                                              % damage status
            uebclogical = logical(uebcvalues);                              % turns into logical array -> 0:fixed nodes 1:moving nodes
            fracture = false;
            tic

            % Staggered Algorithm -----------------------------------------
            istep = 1;
            while istep <= ntsteps
                % automated displacement increment
                if damage_switch ~= 1                                       % elastic region
                    iduinc = duinc;
                elseif damage_switch == 1 && fracture == true               % after fracture
                    iduinc = iduinc*1;
                    if iduinc > duincMAX
                        iduinc = duincMAX;
                    end
                elseif damage_switch == 1 &&  iduinc >= duincMIN            % before fracture
                    iduinc = iduinc/2;
                    if iduinc < duincMIN
                        iduinc = duincMIN;
                        fracture = true;
                    end
                end
                disp(' ');
                dbc = Disp(istep) + iduinc;
                fprintf('Displacement at timestep %d = %1.3e; dinc = %1.3e\n',istep,dbc,iduinc);
                Disp(istep+1) = dbc;
                uebcvalues(uebclogical) = dbc;
                % Time Step Starts
                phievol(:,istep+1) = phievol(:,istep);
                countphi = 0;                                               % phi iteration count
                while true
                    % reiteration algorithm for automated displacement increment
                    if max(phievol(:,istep+1)) > 0.25 && iduinc ~= duincMIN && fracture == false
                        damage_switch = 1;
                        reiter = true;
                        break;
                    end
                    reiter = false;
                    
                    % Elastic Problem -------------------------------------
                    uvdispevol_dummy = uvdispevol(:,istep+1);
                    [uvdispevol(:,istep+1),ReactionF(:,istep+1)] = Solver.elasticPFM(eType,nodes,elements,uebcnodes,uebcvalues,reacnodes,phievol(:,istep+1),YoungsMod,RatPoisson);
                    [strain(:,:,istep+1),Psi_plus] = Solver.Psi_pm(eType,elements,nodes,uvdispevol(:,istep+1),YoungsMod,RatPoisson);
                    Hplus(Hplus<Psi_plus) = Psi_plus(Hplus<Psi_plus);
                    % uvdispevol_error = norm((uvdispevol(:,istep+1)-uvdispevol_dummy)/uvdispevol_dummy);
                    
                    % Damage Problem --------------------------------------
                    phi_dummy = phievol(:,istep+1);
                    [phievol(:,istep+1)] = Solver.damagePFM(eType,nodes,elements,phiebcnodes,phiebcvalues,Hplus,EnergyGc,Charlen);
                    countphi = countphi+1;
                    disp(['Max phi Value: ', num2str(max(phievol(:,istep+1)))]);
                    phistatus = ['Time step : phi-Iteration = ',num2str(istep),':',num2str(countphi)];
                    disp(phistatus);
                    phievol_error = (phievol(:,istep+1)-phi_dummy);
                    % phievol_error2 = norm((phievol(:,istep+1)-phi_dummy)/phi_dummy);

                    if abs(phievol_error) <= tolerance
                        break;
                    end
                end
                % ------------------------------- Update the Nodal Solution
                simdata.uvdispevol = uvdispevol;
                simdata.phievol = phievol;
                simdata.strain = strain;
                simdata.ReactionF = ReactionF;
                simdata.Disp = Disp;
                if damage_switch == 1 && reiter == true
                    continue;
                end
                simdata.istep = istep;
                istep = istep+1;
                elapsedtime = datevec(toc./(60*60*24));
                fprintf("Elapsed time: %d days %d hours %d min %d sec\n",elapsedtime(3),elapsedtime(4),elapsedtime(5),round(elapsedtime(6)))
            end
        end

        function [uvdisp,reaction] = elasticPFM(eType,nodes,elements,uebcnodes,uebcvalues,reacnodes,phi,E,nu)
            %
            %
            %
            if eType == "tri3"
                gw = [0.5 1];
                gp = [1/3 1/3];
            elseif eType == "quad4"
                Ngp = 2;                                                   % No of gauss points in each dimension -- Numerical Integration
                [gp, gw] = Solver.gauss_pw2D(Ngp);                         % Gauss points and weights for 2D numerical integration
            end
            krs = 1e-9;                                                    % Small Residual Stiffness
            Cpsn = (E/(1+nu)/(1-2*nu))*[1-nu nu 0;
                nu 1-nu 0;
                0 0 (1-2*nu)/2];                                           % Plane strain
            % mesh details ------------------------------------------------
            noE = size(elements,1);                                        % number of elements
            noN = size(nodes,1);                                           % number of nodes
            ndof = 2;                                                      % number of dof per node
            % initialize the global matrices ------------------------------
            kmat = zeros(ndof*noN,ndof*noN);
            fvec = zeros(ndof*noN,1);
            % loop over the elements: K-matrix and Body Force -------------
            for iel = 1:noE

                % get the element connectivity
                elcon = elements(iel,:);

                % get nodal coordinates and nodal values
                nds = nodes(elcon,:) ;
                elphi = phi(elcon) ;

                % global index
                gind = zeros(length(elcon)*ndof,1);
                gind(1:2:end) = 2.*elcon-1;
                gind(2:2:end) = 2.*elcon;

                % loop over the gauss points
                for igp = 1:size(gp,1)
                    % get gauss point and respective weight
                    gpt = gp(igp,:);
                    w_igp = gw(igp,1)*gw(igp,2);

                    % get the shape functions and its derivatives
                    [n,dns,dnt] = Solver.getShape2d(gpt,eType);

                    % to locate the point in the physical space
                    % xglobal = n*nds(:,1) ;
                    % yglobal = n*nds(:,2) ;
                    phiglobal = n*elphi;

                    % find the jacobian
                    jac = [dns;dnt]*nds ;

                    % find the derivatives in the physical space
                    dnxy = jac\[dns;dnt];
                    if eType == "tri3"
                        B = [dnxy(1,1) 0 dnxy(1,2) 0 dnxy(1,3) 0;
                            0 dnxy(2,1) 0 dnxy(2,2) 0 dnxy(2,3);
                            dnxy(2,1) dnxy(1,1) dnxy(2,2) dnxy(1,2) dnxy(2,3) dnxy(1,3)];
                    elseif eType == "quad4"
                        B = [dnxy(1,1) 0 dnxy(1,2) 0 dnxy(1,3) 0 dnxy(1,4) 0;
                            0 dnxy(2,1) 0 dnxy(2,2) 0 dnxy(2,3) 0 dnxy(2,4);
                            dnxy(2,1) dnxy(1,1) dnxy(2,2) dnxy(1,2) dnxy(2,3) dnxy(1,3) dnxy(2,4) dnxy(1,4)];
                    end

                    % elemental bilinear and linear form
                    kmat(gind,gind) = kmat(gind,gind) + B'*((1-phiglobal)^2+krs)*Cpsn*B*w_igp*det(jac);

                    % linear form
                    body = [0 0]';
                    if eType == "tri3"
                        Nfun = [n(1) 0 n(2) 0 n(3) 0;
                            0 n(1) 0 n(2) 0 n(3)];
                    elseif eType == "quad4"
                        Nfun = [n(1) 0 n(2) 0 n(3) 0 n(4) 0;
                            0 n(1) 0 n(2) 0 n(3) 0 n(4)];
                    end
                    fvec(gind,1) = fvec(gind,1) + Nfun'*body*w_igp*det(jac);
                end
            end

            % Applying Essential Boundary Conditions ----------------------
            TotalDof = 1:ndof*noN;
            activeDof = setdiff(TotalDof,uebcnodes);

            temp = zeros(size(TotalDof,1),1) ;
            for in=1:length(uebcnodes)
                temp = temp + kmat(:,uebcnodes(in))*uebcvalues(in) ;
            end

            % modified right hand side
            fmod = fvec- temp ;
            fmod = fmod(activeDof,1);
            % modified bilinear form
            kmod = kmat(activeDof,activeDof);
            % solve
            uvactive = kmod\fmod;
            % has to be supplemented with the boundary conditions
            uvdisp(uebcnodes,1) = uebcvalues;
            uvdisp(activeDof,1) = uvactive ;

            % Calculating reactions ---------------------------------------
            reaction = kmat(reacnodes,:)*uvdisp;
        end

        function [phi] = damagePFM(eType,nodes,elements,phiebcnodes,phiebcvalues,Hplus,EnergyGc,Charlen)
            %
            %
            %
            if eType == "tri3"
                gw = [0.5 1];
                gp = [1/3 1/3];
            elseif eType == "quad4"
                Ngp = 2;                                                       % No of gauss points in each dimension -- Numerical Integration
                [gp, gw] = Solver.gauss_pw2D(Ngp);                             % Gauss points and weights for 2D numerical integration
            end
            % mesh details ------------------------------------------------
            noE = size(elements,1);                                        % number of elements
            noN = size(nodes,1);                                           % number of nodes
            ndof = 1;                                                      % number of dof per node
            % initialize the global matrices ------------------------------
            kmat = zeros(ndof*noN,ndof*noN);
            fvec = zeros(ndof*noN,1);
            % loop over the elements: K-matrix and Body Force -------------
            for iel = 1:noE

                % get the element connectivity
                elcon = elements(iel,:);

                % get nodal coordinates and nodal values
                nds = nodes(elcon,:) ;
                % get element values
                elHplus  = Hplus(iel);

                % global index
                gind = elcon;

                % loop over the gauss points
                for igp = 1:size(gp,1)
                    % get gauss point and respective weight
                    gpt = gp(igp,:);
                    w_igp = gw(igp,1)*gw(igp,2);

                    % get the shape functions and its derivatives
                    [n,dns,dnt] = Solver.getShape2d(gpt,eType);

                    % to locate the point in the physical space
                    % xglobal = n*nds(:,1) ;
                    % yglobal = n*nds(:,2) ;

                    % find the jacobian
                    jac = [dns;dnt]*nds ;

                    % find the derivatives in the physical space
                    dnxy = jac\[dns;dnt];

                    % elemental bilinear and linear form
                    kmat(gind,gind) = kmat(gind,gind) + ((dnxy'*EnergyGc*Charlen*dnxy)+(n'*(EnergyGc/Charlen + 2*elHplus)*n))*w_igp*det(jac);

                    % linear form
                    fvec(gind,1) = fvec(gind,1) + n'*2*elHplus*w_igp*det(jac);
                end
            end
            % Applying Essential Boundary Conditions ----------------------
            TotalDof = 1:ndof*noN;
            activeDof = setdiff(TotalDof,phiebcnodes);

            temp = zeros(size(TotalDof,1),1) ;
            for in=1:length(phiebcnodes)
                temp = temp + kmat(:,phiebcnodes(in))*phiebcvalues(in) ;
            end

            %modified right hand side
            fmod = fvec- temp ;
            fmod = fmod(activeDof,1);

            % modified bilinear form
            kmod = kmat(activeDof,activeDof);

            % solve
            phiactive = kmod\fmod;

            % has to be supplemented with the boundary conditions
            phi(phiebcnodes,1) = phiebcvalues;
            phi(activeDof,1) = phiactive ;
            phi(phi<0) = 0;
            phi(phi>1) = 1;
        end

        function [istrain,Psi_plus] = Psi_pm(eType,elements,nodes,uvdisp,E,nu)
            lambda = E*nu/(1+nu)/(1-2*nu);
            mu = E/2/(1+nu);
            if eType == "tri3"
                % gw = [0.5 1];
                gp = [1/3 1/3];
            elseif eType == "quad4"
                Ngp = 1;                                                   % No of gauss points in each dimension -- Numerical Integration
                [gp, gw] = Solver.gauss_pw2D(Ngp);                         % Gauss points and weights for 2D numerical integration
            end
            % mesh details ------------------------------------------------
            noE = size(elements,1);                                        % number of elements
            ndof = 2;                                                      % number of dof per node
            Psi_plus = zeros(noE,1);
            istrain = zeros(3,noE);
            % loop over the elements --------------------------------------
            for iel = 1:noE

                % get the element connectivity
                elcon = elements(iel,:);

                % get nodal coordinates and nodal displacements
                nds = nodes(elcon,:) ;

                % global index
                gind = zeros(length(elcon)*ndof,1);
                gind(1:2:end) = (2.*elcon-1);
                gind(2:2:end) = (2.*elcon);
                eluvdisp = uvdisp(gind) ;

                igp = 1;
                % get gauss point and respective weight
                gpt = gp(igp,:);

                % get the shape functions and its derivatives
                [n,dns,dnt] = Solver.getShape2d(gpt,eType);

                % find the jacobian
                jac = [dns;dnt]*nds ;

                % find the derivatives in the physical space
                dnxy = jac\[dns;dnt];
                if eType == "tri3"
                    B = [dnxy(1,1) 0 dnxy(1,2) 0 dnxy(1,3) 0;
                        0 dnxy(2,1) 0 dnxy(2,2) 0 dnxy(2,3);
                        dnxy(2,1) dnxy(1,1) dnxy(2,2) dnxy(1,2) dnxy(2,3) dnxy(1,3)];
                elseif eType == "quad4"
                    B = [dnxy(1,1) 0 dnxy(1,2) 0 dnxy(1,3) 0 dnxy(1,4) 0;
                        0 dnxy(2,1) 0 dnxy(2,2) 0 dnxy(2,3) 0 dnxy(2,4);
                        dnxy(2,1) dnxy(1,1) dnxy(2,2) dnxy(1,2) dnxy(2,3) dnxy(1,3) dnxy(2,4) dnxy(1,4)];
                end
                strain2d = B*eluvdisp;
                istrain(:,iel) = strain2d;
                strain2dT = [strain2d(1) strain2d(3); strain2d(3) strain2d(2)];
                Psi_plus(iel) = 0.5*lambda*ramp(trace(strain2dT))^2 + mu*trace(rampStrain(strain2dT)^2);
            end
            function [aplus,aminus] = ramp(a)
                aplus = (a+abs(a))/2;
                aminus = (a-abs(a))/2;
            end
            function [Tplus,Tminus] = rampStrain(T)
                [eigdir,eigval] = eig(T);
                Tplus = zeros(size(T));
                Tminus = zeros(size(T));
                for i = 1:length(T)
                    [eigplus,eigminus] = ramp(eigval(i,i));
                    Tplus = Tplus + eigplus*eigdir(:,i)*eigdir(:,i)';
                    Tminus = Tminus + eigminus*eigdir(:,i)*eigdir(:,i)';
                end
            end
        end
        
        function [n,dns,dnt] = getShape2d(pt,eType)
            s = pt(1);
            t = pt(2);
            if eType == "tri3"
                %pt - gauss point
                n = [1-s-t, s, t];
                % derivative is wrt to pt (parameteric coordiante)
                dns = [-1, 1, 0];
                dnt = [-1, 0, 1];
            elseif eType == "quad4"
                %pt - gauss point
                n = [0.25*(1-s)*(1-t), 0.25*(1+s)*(1-t), 0.25*(1+s)*(1+t), 0.25*(1-s)*(1+t)];
                % derivative is wrt to pt (parameteric coordinate)
                dns = [-0.25*(1-t), 0.25*(1-t), 0.25*(1+t), -0.25*(1+t)];
                dnt = [-0.25*(1-s), -0.25*(1+s), 0.25*(1+s), 0.25*(1-s)];
            elseif eType == "quad8"
                %pt - gauss point
                n = [0.25*(1-s)*(t-1)*(1+s+t), 0.25*(1+s)*(t-1)*(t-s+1), 0.25*(1+s)*(1+t)*(s+t-1), 0.25*(s-1)*(t+1)*(s-t+1), 0.5*(1-s^2)*(1-t), 0.5*(1+s)*(1-t^2), 0.5*(1-s^2)*(1+t), 0.5*(1-s)*(1-t^2)];
                % derivative is wrt to pt (parameteric coordinate)
                dns = [-0.25*(t-1)*(2*s+t), -0.25*(t-1)*(2*s-t), 0.25*(t+1)*(2*s+t), 0.25*(t+1)*(2*s-t), 0.5*(-2*s)*(1-t), 0.5*(1-t^2), 0.5*(-2*s)*(1+t), -0.5*(1-t^2)];
                dnt = [0.25*(1-s)*(2*t+s), 0.25*(1+s)*(2*t-s), 0.25*(1+s)*(2*t+s), -0.25*(s-1)*(2*t-s), -0.5*(1-s^2), -0.5*(1+s)*(2*t), 0.5*(1-s^2), -0.5*(1-s)*(2*t)];
            end
        end
        
        % Gauss points and weights ----------------------------------------
        function [gp, gw] = gauss_pw2D(N)
            % HELP - gauss_pw2D (Gauss Points for 2D analysis)
            %   -- generates gauss points and weights for 2D numerical integration
            %   -- input 'N' gives, 'N^2' no of points and weights
            %   -- N: even numbers only [No of gauss points in each dimension]

            % Gauss points 2D
            [x,wx]=Solver.lgwt(N,-1,1);
            [y,wy]=Solver.lgwt(N,-1,1);
            gp = zeros(N^2, 2);
            gw = zeros(N^2, 2);

            row = 1;
            for i = x'
                for j = y'
                    gp(row,:) = [i j];
                    row = row+1;
                end
            end

            row = 1;
            for i = wx'
                for j = wy'
                    gw(row,:) = [i j];
                    row = row+1;
                end
            end
        end
        
        function [x,w] = lgwt(N,a,b)
            % Help - lgwt (Legendre Gauss Weights)
            % -- This script is for computing definite integrals using Legendre-Gauss Quadrature. Computes the Legendre-Gauss nodes and weights on an interval [a,b] with truncation order N
            % -- Suppose you have a continuous function f(x) which is defined on [a,b] which you can evaluate at any x in [a,b]. Simply evaluate it at all of the values contained in the x vector to obtain a vector f. Then compute the definite integral using sum(f.*w);
            % -- Written by Greg von Winckel - 02/25/2004

            % Gauss Points
            N=N-1;
            N1=N+1; N2=N+2;
            xu=linspace(-1,1,N1)';

            % Initial guess
            y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

            % Legendre-Gauss Vandermonde Matrix
            L=zeros(N1,N2);

            % Derivative of LGVM
            Lp=zeros(N1,N2);

            % Compute the zeros of the N+1 Legendre Polynomial
            % using the recursion relation and the Newton-Raphson method
            y0=2;

            % Iterate until new points are uniformly within epsilon of old points
            while max(abs(y-y0))>eps
                L(:,1)=1;
                Lp(:,1)=0;
                L(:,2)=y;
                Lp(:,2)=1;
                for k=2:N1
                    L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
                end
                Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
                y0=y;
                y=y0-L(:,N2)./Lp;
            end

            % Linear map from[-1,1] to [a,b]
            x=(a*(1-y)+b*(1+y))/2;

            % Compute the weights
            w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
        end
    end
end

