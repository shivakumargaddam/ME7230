clear all
close all
clc

% nodal coordiantes
node = [0,0; 1 0; 1,1; 0,1] ;
element = [1 2 3 4];


eType = 4;
nu = 0.499;
E = 1;
% C = eye(3); % Identity
% C = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 1-nu]; % Plane stress
C = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 1-2*nu]; % Plane strain


% generate the elements
numelem = size(element,1); % number of elements
numnode = size(node,1);  % number of nodes per element
ndof = 2; % number of dof per node


% initialize the global matrices
kmat = zeros(ndof*numnode,ndof*numnode);
fvec = zeros(ndof*numnode,1);

% gauss points
order = 2;
if eType == 3
    if order == 1
        w = 0.5;
        pt = [1/3 1/3];
    elseif order == 2
        
    elseif order == 3
        
    end
elseif eType == 4
    if order == 1
        w = [4];
        pt = [0 0];
    elseif order == 2
        w = [1 1 1 1];
        pt = [1/sqrt(3) 1/sqrt(3); -1/sqrt(3) 1/sqrt(3); -1/sqrt(3) -1/sqrt(3); 1/sqrt(3) -1/sqrt(3)];
    elseif order == 3
        
    end
end

% loop over the elements
for iel = 1:numelem
    
    % get the element connectivity
    econ = element(iel,:) ;
    
    % get its coordinates
    nds = node(econ,:) ;
    
    % global index
    gind = [];
    for i = 1:eType
        gind = [gind 2*econ(i)-1 2*econ(i)];
    end
    
    ele = zeros(eType*ndof,eType*ndof);
    % loop over the gauss points
    for igp = 1:size(w,2)
        
        % get the gauss point
        gpt = pt(igp,:);
        
        % get the shape functions and its derivatives
        [n,dns,dnt] = getShape2d(gpt,eType);
        
        % to locate the point in the physical space
        xyglobal = n*nds ;
        
%         % plot the gauss point in the physical space
%         plot(xyglobal(1),xyglobal(2),'bs');
%         hold on;
        
        % find the jacobian
        jac = [dns;dnt]*nds ;
        
        % find the derivatives in the physical space
        dnxy = [1 0 0 0; 0 0 0 1; 0 1 1 0]*[inv(jac)*[dns(1) 0 dns(2) 0 dns(3) 0 dns(4) 0; dnt(1) 0 dnt(2) 0 dnt(3) 0 dnt(4) 0]; inv(jac)*[0 dns(1) 0 dns(2) 0 dns(3) 0 dns(4); 0 dnt(1) 0 dnt(2) 0 dnt(3) 0 dnt(4)]] ;
        
        % elemental bilinear and linear form
        kmat(gind,gind) = kmat(gind,gind) + dnxy'*C*dnxy*w(igp)*det(jac);
        
        % linear form
        b = [0 0]';
        fvec(gind,1) = fvec(gind,1) + [n(1) 0 n(2) 0 n(3) 0 n(4) 0; 0 n(1) 0 n(2) 0 n(3) 0 n(4)]'*b*w(igp)*det(jac);
        
    end
end
kmat
[V,D] = eig(kmat)

for vi = 1:8
    node_ev = node + [V(1,vi),V(2,vi); V(3,vi) V(4,vi); V(5,vi),V(6,vi); V(7,vi),V(8,vi)] ;
        
    figure; hold on;
    set(gca,'FontSize',13)
    plot([node(:,1);node(1,1)],[node(:,2);node(1,2)],'LineWidth',1.5,'DisplayName','undeformed')
    plot([node_ev(:,1);node_ev(1,1)],[node_ev(:,2);node_ev(1,2)],'LineWidth',1.5,'DisplayName','deformed');
    legend('Location','northeast')
    box on
    daspect([1 1 1]);
    set(gcf,'units','pixels','position',[100 100 600 600]);
    xlim([-1,2])
    ylim([-1,2])
    title(['eig-',num2str(vi)])
end

