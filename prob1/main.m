clear all
close all
clc

xo = 0;
x1 = 3;

numpts = 4;
% nodal coordiantes
node = linspace(xo,x1,numpts) ;

% number of nodes
numnode = size(node,2) ;

% generate the elements
numelem = numnode - 1; % number of elements
nnode = 2;  % number of nodes per element
ndof = 1; % number of dof per node

element = zeros(numelem,nnode);
for iel = 1:numelem
    element(iel,:) = [iel iel+1];
end

% % plot the mesh
% figure; hold on;
% plot(node,zeros(numnode,1),'r-o')
% for in = 1:numnode
%     x = node(in);
%     text(x,0.1,num2str(in));
% end
% 
% for iel = 1:numelem
%     econ = element(iel,:) ;
%     nds = node(econ);
%     text(mean(nds),0.01,num2str(iel));
% end

% initialize the global matrices
kmat = zeros(ndof*numnode,ndof*numnode);
fvec = zeros(ndof*numnode,1);

% gauss points
order = 3;
if order == 1
    w = 2;
    pt = 0;
elseif order == 2
    w = [1,1];
    pt = [-1/sqrt(3), 1/sqrt(3)];
elseif order == 3
    w = [5/9 5/9 8/9];
    pt = [-sqrt(3/5) sqrt(3/5) 0];
end

% loop over the elements
for iel = 1:numelem
    
    % get the element connectivity
    econ = element(iel,:) ;
    
    % get its coordinates
    nds = node(econ) ;
    
    % global index
    gind = econ ;
    
    l = 0 ;
    ele = zeros(2,2);
    % loop over the gauss points
    for igp = 1:size(w,2)
        
        % get the gauss point
        gpt = pt(igp);
        
        % get the shape functions and its derivatives
        [n,dn] = getShape1d(gpt);
        
        % to locate the point in the physical space
        xglobal = n*nds' ;
        
        % plot the gauss point in the physical space
        %plot(xglobal(1),0,'bs');
        
        % find the jacobian
        jac = dn*nds' ;
        
        % find the derivatives in the physical space
        dnx = inv(jac)*dn ;
        
        % elemental bilinear and linear form
        kmat(gind,gind) = kmat(gind,gind) + dnx'*dnx*w(igp)*det(jac) ;
        
        % linear form
        b = -4*xglobal ;
        fvec(gind,1) = fvec(gind,1) + n'*b*w(igp)*det(jac);
        
    end
end

% specify the location of boundary conditions
bcdof = [1 numnode];
bcval = [0 52];

TotalDof = 1:ndof*numnode ;
activeDof = setdiff(TotalDof,bcdof);

temp = kmat(:,size(kmat,1))*52 ;

%modified right hand side
fmod = fvec- temp ;

% modified bilinear form
kmod = kmat(activeDof,activeDof) ;

% solve
u = kmod\fmod(activeDof,1)

% has to be supplemented with the boundary conditions
udisp = zeros(ndof*numnode,1) ;
udisp(1) = 0;
udisp(activeDof,1) = u ;
udisp(end) = 52 ;

figure; hold on;
plot(node,udisp,'rs-')

% % plot the analytical solution
% x = linspace(xo,x1,100);
% for i=1:length(x)
%     uana(i,1) = 2/3*x(i)^3 + 34/3*x(i);
% end
% 
% plot(x,uana,'b-')

% loop over the elements for error
L2 = 0;
H1 = 0;
for iel = 1:numelem
    
    % get the element connectivity
    econ = element(iel,:) ;
    
    % get its coordinates
    nds = node(econ) ;
    
    % global index
    gind = econ ;
    % elemental solution
    uele = udisp(gind);
    
    l = 0 ;
    ele = zeros(2,2);
    % loop over the gauss points
    for igp = 1:size(w,2)
        
        % get the gauss point
        gpt = pt(igp);
        
        % get the shape functions and its derivatives
        [n,dn] = getShape1d(gpt);
        
        % to locate the point in the physical space
        xglobal = n*nds' ;
        
        % plot the gauss point in the physical space
        %plot(xglobal(1),0,'bs');
        
        % find the jacobian
        jac = dn*nds' ;
        
        % find the derivatives in the physical space
        dnx = inv(jac)*dn ;
        
        % error for u
        uerror = 2/3*xglobal^3 + 34/3*xglobal - n*uele;
        L2 = L2 + uerror*uerror*det(jac)*w(igp);
        duerror = 2*xglobal^2 + 34/3 - dnx*uele;
        H1 = H1 + duerror*duerror*det(jac)*w(igp);

        % error for u'
        
    end
end

L2 = sqrt(L2)
H1 = sqrt(H1)
h = node(2) - node(1)
