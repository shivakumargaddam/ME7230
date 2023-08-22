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

meshsize = 10;


[node1,element1,gpglobal1,udisp1,dudx1,fvec1] = onedBar (meshsize,1);
[node2,element2,gpglobal2,udisp2,dudx2,fvec2] = onedBar (meshsize,2);
numnode1 = size(node1,2) ;
numnode2 = size(node2,2) ;
numelem1 = size(element1,1) ;
numelem2 = size(element2,1) ;

%% u
figure; hold on;
set(gca,'FontSize',12)
box on
xlim([xstart,xend])
ylim([0 6e-3])
pbaspect([1 1 1]);
xlabel('x')
ylabel('displacement [unit]')
set(gcf,'units','pixels','position',[100 100 600 600]);
title("1D bar: convergence study - 'u'")
plot(node2,zeros(numnode2,1),'.b-','MarkerSize',19,'LineWidth',2,'DisplayName','mesh-qudratic')
plot(node1,zeros(numnode1,1),'.r-','MarkerSize',19,'LineWidth',2,'DisplayName','mesh-linear')
x = linspace(xstart,xend,1000);
for i=1:length(x)
    if x(i) <= x_b
        uana(i,1) = (-1/320)*x(i)*(16*x(i)^4*(k1+k2)-31*k1-k2)/(k1*(k1+k2));
    else
        uana(i,1) = (-1/320)*(16*x(i)^5*(k1+k2)-31*k1*x(i)-k2*x(i)+15*(k1-k2))/(k2*(k1+k2));
    end
end
plot(x,uana,"Color","green",'MarkerSize',7,'LineWidth',7,'DisplayName','Analytical')
plot(node1,udisp1,'-*',"Color","#0072BD",'MarkerSize',9,'LineWidth',1.5,'DisplayName','linear')

udispquad = zeros(size(x));
for ix = 1:length(x)
    for iel = 1:numelem2
        econ = element2(iel,:);
        nds = node2(econ) ;
        if x(ix)>=min(nds) && x(ix)<=max(nds)
            xi = (2*x(ix)-nds(1)-nds(3))/(nds(3)-nds(1));
            udispquad(ix) = 1/2*[-xi+xi^2 2-2*xi^2 xi+xi^2]*udisp2(econ);
        end
    end
end
plot(x,udispquad,"Color","#D95319",'MarkerSize',9,'LineWidth',1.5,'DisplayName','quadratic')
plot(node2,udisp2,'*',"Color","#D95319")
legend('Location','northeast')

%% du/dx
figure; hold on;
set(gca,'FontSize',12)
box on
xlim([xstart,xend])
% ylim([0 6e-3])
pbaspect([1 1 1]);
xlabel('x')
ylabel('du/dx')
set(gcf,'units','pixels','position',[100 100 600 600]);
title("1D bar: convergence study - 'du/dx'")
plot(node2,zeros(numnode2,1),'.b-','MarkerSize',19,'LineWidth',1.5,'DisplayName','mesh-qudratic')
plot(node1,zeros(numnode1,1),'.r-','MarkerSize',19,'LineWidth',1.5,'DisplayName','mesh-linear')
x = linspace(xstart,xend,1000);
for i=1:length(x)
    if x(i) <= x_b
        duana(i,1) = (-1/320/k1/(k1+k2))*(80*x(i)^4*(k1+k2)-31*k1-k2);
    else
        duana(i,1) = (-1/320/k2/(k1+k2))*(80*x(i)^4*(k1+k2)-31*k1-k2);
    end
end
plot(x,duana,"Color","green",'MarkerSize',7,'LineWidth',7,'DisplayName','Analytical')
dudxlin = zeros(size(x));
for ix = 1:length(x)
    for iel = 1:numelem1
        econ = element1(iel,:);
        nds = node1(econ) ;
        if x(ix)>=min(nds) && x(ix)<=max(nds)
            dudxlin(ix) = dudx1(3*iel);
        end
    end
end
plot(x,dudxlin,'.',"Color","#0072BD",'MarkerSize',5,'LineWidth',1.5,'DisplayName','linear')
plot(gpglobal2,dudx2,"Color","#D95319",'MarkerSize',5,'LineWidth',1.5,'DisplayName','quadratic')
% plot(node2,udisp2,'*',"Color","#D95319")
legend('Location','northeast')

%% force vector
figure; hold on;
set(gca,'FontSize',12)
box on
xlim([xstart,xend])
ylim([0 Inf])
pbaspect([1 1 1]);
xlabel('x')
ylabel('q(x)/force vector')
legend('Location','northwest')
set(gcf,'units','pixels','position',[100 100 600 600]);
title("1D bar: convergence study - force vector")
plot(node2,zeros(numnode2,1),'.b-','MarkerSize',19,'LineWidth',2,'DisplayName','mesh-qudratic')
plot(node1,zeros(numnode1,1),'.r-','MarkerSize',19,'LineWidth',2,'DisplayName','mesh-linear')
x = linspace(xstart,xend,1000);
plot(x,x.^3,"Color","green",'MarkerSize',7,'LineWidth',7,'DisplayName','Analytical')
plot(node1,fvec1,"Color","#0072BD",'MarkerSize',9,'LineWidth',1.5,'DisplayName','linear')
plot(node2,fvec2,"Color","#D95319",'MarkerSize',9,'LineWidth',1.5,'DisplayName','quadratic')




function [node,element,gpglobal,udisp,dudx,fvec] = onedBar (divisions,order)
% Input
k1 = 1;
k2 = 10;
Length = 1;
x_b = 0.5;


% Variables
xstart = 0;
xend = Length;

% Variables - accuracy of solution
numpts = divisions + 1; % no of divisions = numpts-1
eleOrder = order;
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


%% 

% gauss points
if nofgp == 1
    w = 2;
    pt = 0;
elseif nofgp == 2
    w = [1,1];
    pt = [-1/sqrt(3), 1/sqrt(3)];
elseif nofgp == 3
    w = [5/9 8/9 5/9];
    pt = [-sqrt(3/5) 0 sqrt(3/5)];
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
        kmat(gind,gind) = kmat(gind,gind) + k*dnx'*dnx*w(igp)*det(jac) ;
        
        % linear form
        b = xglobal^3 ;
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
fmod = fvec- temp 

% modified bilinear form
kmod = kmat(activeDof,activeDof) ;

% solve
u = kmod\fmod(activeDof,1);

% has to be supplemented with the boundary conditions
udisp(ebcdof) = ebcval;
udisp(activeDof,1) = u ;

%% plot the analytical solution
if false
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

    x = linspace(xstart,xend,1000);
    for i=1:length(x)
        if x(i) <= x_b
            uana(i,1) = (-1/320)*x(i)*(16*x(i)^4*(k1+k2)-31*k1-k2)/(k1*(k1+k2));
        else
            uana(i,1) = (-1/320)*(16*x(i)^5*(k1+k2)-31*k1*x(i)-k2*x(i)+15*(k1-k2))/(k2*(k1+k2));
        end
    end
    plot(x,uana,"Color","#77AC30",'MarkerSize',7,'LineWidth',3,'DisplayName','Analytical')
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


%% dudx estimation
dudx = [];
gpglobal = [];

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
        gpglobal = [gpglobal xglobal];

        % find the jacobian
        jac = dn*nds' ;
        
        % find the derivatives in the physical space
        dnx = inv(jac)*dn ;
        
        dudx = [dudx dnx*uele];
    end
end

% length(dudx) == numelem*nofgp
% 
% length(dudx) == length(gpglobal)

end

