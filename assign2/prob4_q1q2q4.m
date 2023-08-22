clear all
close all
clc

%% 1
FGPswitch = 0;
k = 9;
natfand = (2.*k+1).*pi./2;

meshsize1 = 10:40;
whwa1 = zeros(size(meshsize1));
for imeshsize = 1:length(meshsize1)
    [natf] = axialvibrations (meshsize1(imeshsize),1,FGPswitch,0);
    whwa1(imeshsize) = natf(10)/natfand;
end

meshsize2 = 5:40;
whwa2 = zeros(size(meshsize2));
for imeshsize = 1:length(meshsize2)
    [natf] = axialvibrations (meshsize2(imeshsize),2,FGPswitch,0);
    whwa2(imeshsize) = natf(10)/natfand;
end


figure; hold on;
set(gca,'FontSize',12)
box on
% xlim([0,1])
% ylim([0 6e-3])
pbaspect([1 1 1]);
xlabel('x')
ylabel('\omega^h/\omega^a')
set(gcf,'units','pixels','position',[100 100 600 600]);
title("Axial vibration of an elastic bar")
plot(meshsize1,whwa1,'.b-','MarkerSize',19,'LineWidth',2,'DisplayName','mesh-linear')
plot(meshsize2,whwa2,'.r-','MarkerSize',19,'LineWidth',2,'DisplayName','mesh-quadratic')
yline(1.05,'-','LineWidth',2,'DisplayName','1.05');
legend

%% 2
FGPswitch = 0;
k = 0:9;
natfand = (2.*k+1).*pi./2;


[natf] = axialvibrations (28,1,FGPswitch,0);
whwak1 = natf(1:10)./natfand';


[natf] = axialvibrations (10,2,FGPswitch,0);
whwak2 = natf(1:10)./natfand';

figure; hold on;
set(gca,'FontSize',12)
box on
% xlim([0,1])
% ylim([0 6e-3])
pbaspect([1 1 1]);
xlabel('mode number (k)')
ylabel('\omega^h/\omega^a')
set(gcf,'units','pixels','position',[100 100 600 600]);
title("Axial vibration of an elastic bar: \omega^h/\omega^a vs. k")
plot(k,whwak1,'.b-','MarkerSize',19,'LineWidth',2,'DisplayName','mesh-linear')
plot(k,whwak2,'.r-','MarkerSize',19,'LineWidth',2,'DisplayName','mesh-quadratic')
% yline(1.05,'-','LineWidth',2,'DisplayName','1.05');
legend('Location','northwest')

%% 4
FGPswitch = 1;
n = 0:0.1:1;

fff = zeros(size(n));
for in = 1:length(n)
    [natf] = axialvibrations (100,2,FGPswitch,n(in));
    fff(in) = natf(1);
end

figure; hold on;
set(gca,'FontSize',12)
box on
xlim([0,1])
% ylim([0 6e-3])
pbaspect([1 1 1]);
xlabel('material gradient index (n)')
ylabel('\omega^h - mode 1')
set(gcf,'units','pixels','position',[100 100 600 600]);
title("Axial vibration of an elastic bar: FGM")
plot(n,fff,'.b-','MarkerSize',19,'LineWidth',2,'DisplayName','mesh-linear')


function [natf] = axialvibrations (divisions,order,FGPswitch,imgi)
% Material Variables
A = 1;
Length = 1;
FGP = FGPswitch;
Ec = 348e9;
Em = 201e9;
rhoc = 2370;
rhom = 8166;
nmgi = imgi;

% Variables
xstart = 0;
xend = Length;


% Variables - accuracy of solution
numpts = divisions + 1; % no of divisions = numpts-1
eleOrder = order;
if eleOrder == 1
    nofgp = 3;
elseif eleOrder == 2
    nofgp = 5;
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

% initialize the global matrices
kmat = zeros(ndof*numnode,ndof*numnode);
mmat = zeros(ndof*numnode,ndof*numnode);
fvec = zeros(ndof*numnode,1);

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
        
        if FGP == 0
            E = 1;
            rho = 1;
        elseif FGP == 1
            Vc = (1-(xglobal/Length))^nmgi;
            Vm = 1-Vc;
            E = Ec*Vc + Em*Vm;
            rho = rhoc*Vc + rhom*Vm;
        end
        
        kmat(gind,gind) = kmat(gind,gind) + (dnx'*E*A*dnx)*w(igp)*det(jac) ;
        mmat(gind,gind) = mmat(gind,gind) + n'*rho*A*n*w(igp)*det(jac) ;
        
    end
end
%% Applying Boundary Conditions

% specify the location of essential boundary conditions
ebcdof = [1];

TotalDof = 1:ndof*numnode ;
activeDof = setdiff(TotalDof,ebcdof);

% modified bilinear form
kmod = kmat(activeDof,activeDof) ;
mmod = mmat(activeDof,activeDof) ;

% solve
natf = sqrt(eig(kmod,mmod));
end