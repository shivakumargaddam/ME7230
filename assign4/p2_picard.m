clear all
close all
clc

% Input
Length = 0.4;
Area = 100e-6;
epsilon1 = 0.15;
Eo = 70e9;
E1 = 49e9;
alpha = 1/epsilon1*(1-E1/Eo);
uguess = [0;0;0;0];

% Variables
xstart = 0;
xend = Length;

% Variables - accuracy of solution
numpts = 4; % no of divisions = numpts-1
eleOrder = 1;
if eleOrder == 1
    nofgp = 3;
elseif eleOrder == 2
    nofgp = 3;
end

%% Mesh Generation
% nodal coordinates
if eleOrder == 1
    node = linspace(xstart,xend,numpts) ;
elseif eleOrder == 2
    node = linspace(xstart,xend,2*numpts-1);
end

% number of nodes
numnode = size(node,2) ;

% generate the elements
if eleOrder == 1
    numelem = numnode - 1; % number of elements
    nnode = 2;  % number of nodes per element
    ndof = 1; % number of dof per node
elseif eleOrder == 2
    numelem = (numnode-1)*0.5; % number of elements
    nnode = 3;  % number of nodes per element
    ndof = 1; % number of dof per node
end

element = zeros(numelem,nnode);
if eleOrder == 1
    for iel = 1:numelem
        element(iel,:) = [iel iel+1];
    end
elseif eleOrder == 2
    for iel = 1:numelem
        element(iel,:) = [2*iel-1 2*iel 2*iel+1];
    end
end

% gauss points
if nofgp == 1
    w = 2;
    pt = 0;
elseif nofgp == 2
    w = [1,1];
    pt = [-1/sqrt(3), 1/sqrt(3)];
elseif nofgp == 3
    w = [5/9 5/9 8/9];
    pt = [-sqrt(3/5) sqrt(3/5) 0];
elseif nofgp == 4
    w = [(18-sqrt(30))/36 (18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36];
    pt = [-sqrt(3/7+2/7*sqrt(6/5)) sqrt(3/7+2/7*sqrt(6/5)) -sqrt(3/7-2/7*sqrt(6/5)) sqrt(3/7-2/7*sqrt(6/5))];
elseif nofgp == 5
    w = [(322-13*sqrt(70))/900 (322-13*sqrt(70))/900 (322+13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225];
    pt = [-1/3*sqrt(5+2*sqrt(10/7)) 1/3*sqrt(5+2*sqrt(10/7)) -1/3*sqrt(5-2*sqrt(10/7)) 1/3*sqrt(5-2*sqrt(10/7)) 0];
end

iter = 1;
while true
    % initialize the global matrices
    kmat = zeros(ndof*numnode,ndof*numnode);
    fvec = zeros(ndof*numnode,1);
    udisp = zeros(ndof*numnode,1);
    
    % loop over the elements
    for iel = 1:numelem
        
        % get the element connectivity
        econ = element(iel,:) ;
        
        % get its coordinates
        nds = node(econ) ;
        
        % global index
        gind = econ ;
        
        % loop over the gauss points
        for igp = 1:size(w,2)
            
            % get the gauss point
            gpt = pt(igp);
            
            % get the shape functions and its derivatives
            [n,dn] = getShape1d(gpt,eleOrder);
            
            % to locate the point in the physical space
            xglobal = n*nds' ;
    
            % find the jacobian
            jac = dn*nds' ;
            
            % find the derivatives in the physical space
            dnx = inv(jac)*dn ;

            E_u = Eo*(1-alpha*dnx*uguess(econ));
            
            % elemental bilinear form
            kmat(gind,gind) = kmat(gind,gind) + (dnx'*Area*E_u*dnx)*w(igp)*det(jac) ;
            
            % linear form
            b = 0;
            fvec(gind,1) = fvec(gind,1) + n'*b*w(igp)*det(jac);
        end
    end
    
    %% Applying Boundary Conditions
    % specify the location of natural boundary conditions
    nbcdof = [numnode];
    nbcval = [800e3];
    fvec(nbcdof,1) = fvec(nbcdof,1) + nbcval;
    
    % specify the location of essential boundary conditions
    ebcdof = [1];
    ebcval = [0];
    
    TotalDof = 1:ndof*numnode ;
    activeDof = setdiff(TotalDof,ebcdof);
    
    temp = zeros(size(TotalDof,1),1) ;
    for in=1:length(ebcdof)
        temp = temp + kmat(:,ebcdof(in))*ebcval(in) ;
    end
    
    %modified right hand side
    fmod = fvec- temp ;
    
    % modified bilinear form
    kmod = kmat(activeDof,activeDof) ;
    
    % solve
    u = kmod\fmod(activeDof,1);
    
    % has to be supplemented with the boundary conditions
    udisp(ebcdof) = ebcval;
    udisp(activeDof,1) = u ;
    
    fprintf("iteration: %d\n",iter)
    if abs(udisp-uguess) < 1e-9
        break;
    else
        uguess = udisp;
        disp(udisp*1e3)
        iter = iter +1;
    end
end
udisp*1e3

