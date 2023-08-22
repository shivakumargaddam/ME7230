tic
clear all
close all
clc

% Material Variables
D = 259.2*(1e-4);
alpha = 0.002125;
depth = 195*(1e-2);
c0 = 0.21;

% Variables
xstart = 0;
xend = depth;
ttime = 5;
dtime = .001;
% type of time integration theta = [0(forward eluler) 1(backward euler) 0.5(crank-nicklson)]
theta = 0; 

% Variables - accuracy of solution
numpts = 128; % no of divisions = numpts-1
eleOrder = 1;
nofgp = 2;

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
end


% initialize the global matrices
kmat = zeros(ndof*numnode,ndof*numnode);
mmat = zeros(ndof*numnode,ndof*numnode);
fvec = zeros(ndof*numnode,1);
conc = zeros(ndof*numnode,1+fix(ttime/dtime));

% initial condition
conc(:,1) = c0;

tstep = 0:dtime:ttime;
%% plot the analytical solution

figure; hold on;
x = node;
x(1) = 1e-5;
t = tstep;
for it = round(linspace(1,length(tstep),10))
    for ix=1:length(x)
        cana(ix,it) = c0 - alpha*t(it) + alpha*((t(it)+(x(ix)^2/2/D))*erfc(x(ix)/2/sqrt(D*t(it))) - x(ix)*(t(it)/pi/D)^0.5*exp(-x(ix)^2/4/D/t(it)));
    end
    p1 = plot(x*100,cana(:,it),"-",'Color',"green",'LineWidth',6);
end

p2 = plot(node*100,conc(:,1),'r','LineWidth',1.5,'DisplayName','t=0');
text(node(end)*100,conc(end,1),'\leftarrow t = 0')

% loop over the time

for istep = 2:length(tstep)
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
            
            kmat(gind,gind) = kmat(gind,gind) + D*dnx'*dnx*w(igp)*det(jac) ;
            mmat(gind,gind) = mmat(gind,gind) + n'*n*w(igp)*det(jac) ;
            fvec(gind,1) = fvec(gind,1) - n'*alpha*w(igp)*det(jac);
        end
    end
    

    %% Applying Boundary Conditions
    % specify the location of natural boundary conditions
    nbcdof = [numnode];
    nbcval = [0];
    fvec(nbcdof,1) = fvec(nbcdof,1) - nbcval;

    % specify the location of essential boundary conditions
    ebcdof = [1];
    ebcval = [c0];
    
    TotalDof = 1:ndof*numnode ;
    activeDof = setdiff(TotalDof,ebcdof);
    
    if theta == 0
        A = mmat;
        B = dtime*fvec+mmat*conc(:,istep-1)-dtime*kmat*conc(:,istep-1);
    else
        A = mmat+dtime*theta*kmat;
        B = dtime*fvec+mmat*conc(:,istep-1)-(1-theta)*dtime*kmat*conc(:,istep-1);
    end

    etemp = A(:,ebcdof)*ebcval ;

    %modified right hand side
    Bmod = B - etemp ;

    % modified bilinear form
    Amod = A(activeDof,activeDof) ;

    % solve
    c = Amod\Bmod(activeDof,1);

    % has to be supplemented with the boundary conditions
    conc(ebcdof,istep) = ebcval;
    conc(activeDof,istep) = c ;

    if ismember(istep,round(linspace(1,length(tstep),10)))
        plot(node*100,conc(:,istep),'b','MarkerSize',7,'LineWidth',1.5)
        text(node(end)*100,conc(end,istep),"\leftarrow t = "+num2str(round(tstep(istep),2))+" h.")
    end
end

set(gca,'FontSize',12)
box on
xlim([xstart,xend].*100)
ylim([0.198 0.212])
pbaspect([1 1 1]);
xlabel('depth [cm]')
ylabel('concentration')
set(gcf,'units','pixels','position',[100 100 600 600]);
title("Oxygen diffusion in soil - explicit")
legend([p1 p2],{'Analytical','t=0'})
% legend

if theta == 0
    disp("Forward-Euler(Explicit)")
    savefig("FE")
elseif theta == 0.5
    disp("Cranck-Nicolson(Implicit)")
    savefig("CN")
elseif theta == 1
    disp("Backward-Euler(Implicit)")
    savefig("BE")
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
    cele = conc(gind,end);
    
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
        cerror = c0 - alpha*ttime + alpha*((ttime+(xglobal^2/2/D))*erfc(xglobal/2/sqrt(D*ttime)) - xglobal*(ttime/pi/D)^0.5*exp(-xglobal^2/4/D/ttime)) - n*cele;
        L2 = L2 + cerror*cerror*det(jac)*w(igp);
    end
end

L2 = sqrt(L2)
h = node(2) - node(1);
toc






























