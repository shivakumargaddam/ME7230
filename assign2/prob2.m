clear all
close all
clc

% Input
k1 = 1;
k2 = 10;
Length = 1;
x_b = 0.5;


% Variables
xstart = 0;
xend = Length;

% Variables - accuracy of solution
numpts = 11; % no of divisions = numpts-1
eleOrder = 2;
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

% % plot the mesh
% figure; hold on;
% box on;
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

%% 

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
        
        % elemental bilinear form
        if xglobal < x_b
            k = k1;
        elseif xglobal > x_b
            k = k2;
        end
        kmat(gind,gind) = kmat(gind,gind) + k*(dnx'*dnx)*w(igp)*det(jac) ;
        
        % linear form
        b = xglobal^3 ;
        % b = 0;
        fvec(gind,1) = fvec(gind,1) + n'*b*w(igp)*det(jac);
        
    end
end
%% Applying Boundary Conditions
% specify the location of natural boundary conditions
nbcdof = [];
nbcval = [];

% specify the location of essential boundary conditions
ebcdof = [1 numnode];
ebcval = [0 0];

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

%% plot the analytical solution
if true
    figure; hold on;
    set(gca,'FontSize',12)
    box on
    xlim([min(node),max(node)])
    pbaspect([1 1 1]);
    xlabel('1D-domain [unit]')
    ylabel('displacement [unit]')
    legend('Location','northeast')
    set(gcf,'units','pixels','position',[100 100 600 600]);
    title("One dimensional bar - convergence study")
    plot(node,zeros(numnode,1),'r-o')
    x = linspace(xstart,xend,1000);
    for i=1:length(x)
        if x(i) <= x_b
            uana(i,1) = (-1/320)*x(i)*(16*x(i)^4*(k1+k2)-31*k1-k2)/(k1*(k1+k2));
        else
            uana(i,1) = (-1/320)*(16*x(i)^5*(k1+k2)-31*k1*x(i)-k2*x(i)+15*(k1-k2))/(k2*(k1+k2));
        end
    end
    plot(x,uana,"Color","green",'MarkerSize',7,'LineWidth',9,'DisplayName','Analytical')
    if eleOrder == 1
        plot(node,udisp,"Color","#0072BD",'MarkerSize',7,'LineWidth',2,'DisplayName','FEM solution')
    elseif eleOrder == 2
        udispquad = zeros(size(x));
        for ix = 1:length(x)
            for iel = 1:numelem
                econ = element(iel,:);
                nds = node(econ) ;
                if x(ix)>=min(nds) && x(ix)<=max(nds)
                    xi = (2*x(ix)-nds(1)-nds(3))/(nds(3)-nds(1));
                    udispquad(ix) = 1/2*[-xi+xi^2 2-2*xi^2 xi+xi^2]*udisp(econ);
                end
            end
        end
        plot(x,udispquad,"Color","#0072BD",'MarkerSize',7,'LineWidth',2,'DisplayName','FEM solution')
    end
end

%% L2 and H1 
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
        [n,dn] = getShape1d(gpt,eleOrder);
        
        % to locate the point in the physical space
        xglobal = n*nds' ;

        % find the jacobian
        jac = dn*nds' ;
        
        % find the derivatives in the physical space
        dnx = inv(jac)*dn ;
        
        % error for u
        if xglobal <= x_b
            uerror = (-1/320)*xglobal*(16*xglobal^4*(k1+k2)-31*k1-k2)/(k1*(k1+k2)) - n*uele;
        else
            uerror = (-1/320)*(16*xglobal^5*(k1+k2)-31*k1*xglobal-k2*xglobal+15*(k1-k2))/(k2*(k1+k2)) - n*uele;
        end
        L2 = L2 + uerror*uerror*det(jac)*w(igp);
        if xglobal <= x_b
            duerror =  (-1/320/k1/(k1+k2))*(80*xglobal^4*(k1+k2)-31*k1-k2) - dnx*uele;
        else
            duerror =  (-1/320/k2/(k1+k2))*(80*xglobal^4*(k1+k2)-31*k1-k2)- dnx*uele;
        end
        H1 = H1 + duerror*duerror*det(jac)*w(igp);
    end
end

L2 = sqrt(L2)
H1 = sqrt(H1)
h = node(2) - node(1)

%% Error at nodes and integration points
x = node;
uana = zeros(size(udisp));
for i=1:length(x)
    if x(i) <= x_b
        uana(i,1) = (-1/320)*x(i)*(16*x(i)^4*(k1+k2)-31*k1-k2)/(k1*(k1+k2));
    else
        uana(i,1) = (-1/320)*(16*x(i)^5*(k1+k2)-31*k1*x(i)-k2*x(i)+15*(k1-k2))/(k2*(k1+k2));
    end
end
udisp-uana
unodeerror = norm(udisp-uana)
