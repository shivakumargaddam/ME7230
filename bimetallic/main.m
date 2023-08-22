clear all
% close all
% clc

x0 = 0;
x1 = 1;
k1 = 1;
k2 = 10;

numpts = 40;
oelement = 1;


% nodal coordiantes
if oelement == 1
    node = linspace(x0,x1,numpts) ;
elseif oelement == 2
    node = linspace(x0,x1,2*numpts-1) ;
end

% number of nodes
numnode = size(node,2) ;

% generate the elements
if oelement == 1
    numelem = numnode - 1; % number of elements
    nnode = 2;  % number of nodes per element
    ndof = 1; % number of dof per node
elseif oelement == 2
    numelem = (numnode-1)*0.5; % number of elements
    nnode = 3;  % number of nodes per element
    ndof = 1; % number of dof per node
end

element = zeros(numelem,nnode);
if oelement == 1
    for iel = 1:numelem
        element(iel,:) = [iel iel+1];
    end
elseif oelement == 2
    for iel = 1:numelem
        element(iel,:) = [2*iel-1 2*iel 2*iel+1];
    end
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
order = 2;
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
    
    % loop over the gauss points
    for igp = 1:size(w,2)
        
        % get the gauss point
        gpt = pt(igp);
        
        % get the shape functions and its derivatives
        [n,dn] = getShape1d(gpt,oelement);
        
        % to locate the point in the physical space
        xglobal = n*nds' ;
        
        % plot the gauss point in the physical space
        %plot(xglobal(1),0,'bs');
        
        % find the jacobian
        jac = dn*nds' ;
        
        % find the derivatives in the physical space
        dnx = inv(jac)*dn ;
        
        % elemental bilinear and linear form
        if xglobal < 0.5
            k = k1;
        elseif xglobal > 0.5
            k = k2;
        end
        kmat(gind,gind) = kmat(gind,gind) + k*dnx'*dnx*w(igp)*det(jac) ;

        
        % linear form
        b = xglobal^3 ;
        fvec(gind,1) = fvec(gind,1) + n'*b*w(igp)*det(jac);
        
    end
end

% kmat
% fvec

% specify the location of boundary conditions
bcdof = [1 numnode];
bcval = [0 0];

TotalDof = 1:ndof*numnode ;
activeDof = setdiff(TotalDof,bcdof);

temp = kmat(:,size(kmat,1))*0 ;

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
udisp(end) = bcval(end) ;

% figure; hold on;
set(gca,'FontSize',12)
plot(node,udisp,'bs-','MarkerSize',7,'LineWidth',1)

% plot the analytical solution
x = linspace(x0,x1,100);
for i=1:length(x)
%     uana(i,1) = 2/3*x(i)^3 + 34/3*x(i);
%     uana(i,1) = (1/U)*(x(i)-(1-exp(U*x(i)/kappa))/(1-exp(U/kappa)));
    if x(i) <= 0.5
        uana(i,1) = (-1/320)*x(i)*(16*x(i)^4*(k1+k2)-31*k1-k2)/(k1*(k1+k2));
    else
        uana(i,1) = (-1/320)*(16*x(i)^5*(k1+k2)-31*k1*x(i)-k2*x(i)+15*(k1-k2))/(k2*(k1+k2));
    end
end

plot(x,uana,'b-','MarkerSize',7,'LineWidth',1)
box on
pbaspect([1 1 1]);

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
    
    % loop over the gauss points
    for igp = 1:size(w,2)
        
        % get the gauss point
        gpt = pt(igp);
        
        % get the shape functions and its derivatives
        [n,dn] = getShape1d(gpt,oelement);
        
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
