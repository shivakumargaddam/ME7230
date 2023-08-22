function [L2,H1] = cylipipe (nodes,elements,ubcnode,vbcnode,belements,nu,E,P,method)
% --- Input Arguments ---
% method - 1: normal, 2: selective reduced integration, 3: u-p formulation
% E - young's modulus
% P - pressure
% nu - poisson's ratio
% --- Output Arguments ---
% L2 norm, H1 norm


C = (E/(1+nu)/(1-2*nu))*[1-nu nu 0; nu 1-nu 0; 0 0 (1-2*nu)/2]; % Plane strain

% mesh details
numelem = size(elements,1); % number of elements
numbelem = size(belements,1); % number of boundary elements
numnode = size(nodes,1);  % number of nodes
ndof = 2; % number of dof per node

% initialize the global matrices
kmat = zeros(ndof*numnode,ndof*numnode);
fvec = zeros(ndof*numnode,1);

%% loop over the elements: K-matrix and Body Force
for iel = 1:numelem
    % get the gauss points and weights
    eType = 8;
    order = 3;
    [w,pt] = gaussptwt (eType,order);

    % get the element connectivity
    econ = elements(iel,:);

    % get its coordinates
    nds = nodes(econ,:) ;

    % global index
    gind = [];
    for i = 1:length(econ)
        gind = [gind 2*econ(i)-1 2*econ(i)];
    end

    % loop over the gauss points
    for igp = 1:size(w,2)

        % get the gauss point
        gpt = pt(igp,:);

        % get the shape functions and its derivatives
        [n,dns,dnt] = getShape2d(gpt,eType);

        % to locate the point in the physical space
        xglobal = n*nds(:,1) ;
        yglobal = n*nds(:,2) ;

        % find the jacobian
        jac = [dns;dnt]*nds ;

        % find the derivatives in the physical space
        dnxy = inv(jac)*[dns;dnt];
        B = [dnxy(1,1) 0 dnxy(1,2) 0 dnxy(1,3) 0 dnxy(1,4) 0 dnxy(1,5) 0 dnxy(1,6) 0 dnxy(1,7) 0 dnxy(1,8) 0;
            0 dnxy(2,1) 0 dnxy(2,2) 0 dnxy(2,3) 0 dnxy(2,4) 0 dnxy(2,5) 0 dnxy(2,6) 0 dnxy(2,7) 0 dnxy(2,8) ;
            dnxy(2,1) dnxy(1,1) dnxy(2,2) dnxy(1,2) dnxy(2,3) dnxy(1,3) dnxy(2,4) dnxy(1,4) dnxy(2,5) dnxy(1,5) dnxy(2,6) dnxy(1,6) dnxy(2,7) dnxy(1,7) dnxy(2,8) dnxy(1,8)];

        % elemental bilinear and linear form
        kmat(gind,gind) = kmat(gind,gind) + B'*C*B*w(igp)*det(jac);

        % linear form
        body = [0 0]';
        bforce = [n(1) 0 n(2) 0 n(3) 0 n(4) 0 n(5) 0 n(6) 0 n(7) 0 n(8) 0;
            0 n(1) 0 n(2) 0 n(3) 0 n(4) 0 n(5) 0 n(6) 0 n(7) 0 n(8)]'*body*w(igp)*det(jac);
        fvec(gind,1) = fvec(gind,1) + bforce;
    end
end

%% loop over the elements: Traction Force aka boundary integral
for iel = 1:numbelem
    % get the gauss points and weights
    eType = 1;
    order = 2;
    [w,pt] = gaussptwt (eType,order);

    % get the element connectivity
    econ = belements(iel,:) ;

    % get its coordinates
    enodes = nodes(econ,:);
    elength = norm(enodes(1,:)-enodes(2,:));
    normalv = [((enodes(2,2)-enodes(1,2))/elength); -((enodes(2,1)-enodes(1,1))/elength)];
    nds = [-elength/2 elength/2];

    % global index
    gind = [];
    for i = 1:length(econ)
        gind = [gind 2*econ(i)-1 2*econ(i)];
    end

    % loop over the gauss points
    for igp = 1:size(w,2)
        % get the gauss point
        gpt = pt(igp);
        
        % get the shape functions and its derivatives
        [n,dn] = getShape1d(gpt,eType);
        
        % to locate the point in the physical space
        xglobal = n*nds' ;

        % find the jacobian
        jac = dn*nds' ;
        
        % find the derivatives in the physical space
        % dnx = inv(jac)*dn ;
        
        % elemental bilinear form
        tforce = [n(1) 0 n(2) 0; 0 n(1) 0 n(2)]'*P*eye(2)*normalv*w(igp)*det(jac);
        if iel == 1
            tforce = tforce + [-tforce(1);tforce(2);0;0];
        elseif iel == numbelem
            tforce = tforce + [0;0;tforce(3);-tforce(4)];
        end
        fvec(gind,1) = fvec(gind,1) + tforce;
    end
end
%% Applying Essential Boundary Conditions
% specify the location of essential boundary conditions
uebcdof = ubcnode*2-1;
vebcdof = vbcnode*2;
uebcval = zeros(size(uebcdof));
vebcval = zeros(size(vebcdof));
% u2ebcdof = unique(belements)*2-1;
% v2ebcdof = unique(belements)*2;
% u2ebcval = zeros(size(u2ebcdof))+nodes(unique(belements),1);
% v2ebcval = zeros(size(v2ebcdof))+nodes(unique(belements),2);

TotalDof = 1:ndof*numnode;
activeDof = setdiff(TotalDof,[uebcdof;vebcdof]');

temp = zeros(size(TotalDof,1),1) ;
for in=1:length(uebcdof)
    temp = temp + kmat(:,uebcdof(in))*uebcval(in) ;
end
for in=1:length(vebcdof)
    temp = temp + kmat(:,vebcdof(in))*vebcval(in) ;
end
% for in=1:length(uebcdof)
%     temp = temp + kmat(:,u2ebcdof(in))*u2ebcval(in) ;
% end
% for in=1:length(vebcdof)
%     temp = temp + kmat(:,v2ebcdof(in))*v2ebcval(in) ;
% end

%modified right hand side
fmod = fvec- temp ;

% modified bilinear form
kmod = kmat(activeDof,activeDof);

% solve
u = kmod\fmod(activeDof,1);

% has to be supplemented with the boundary conditions
udisp(uebcdof,1) = uebcval;
udisp(vebcdof,1) = vebcval;
% udisp(u2ebcdof,1) = u2ebcval;
% udisp(v2ebcdof,1) = v2ebcval;
udisp(activeDof,1) = u ;

%% Plot solution
aradius = 0.002;
bradius = 0.001;
if true
    meshfig = figure('Name','Mesh','NumberTitle','off');
    grid off;
    axis on;
    xlabel ('x-axis [mm]');
    ylabel ('y-axis [mm]');
    daspect([1 1 1]);
    set(gca, 'color', 'none');                                              % To remove background
    set(gcf,'units','pixels','position',[500 300 800 600]);                 % To change the size of the figure
    patch('Faces',elements,'Vertices',nodes*1e3,'FaceColor','cyan','FaceAlpha',1);
    hold on
    plot(nodes(vbcnode,1)*1e3,nodes(vbcnode,2)*1e3,'*b');
    plot(nodes(ubcnode,1)*1e3,nodes(ubcnode,2)*1e3,'*r');
    % plot(nodes(fbcnode,1)*1e3,nodes(fbcnode,2)*1e3,'*k');
    patch('Faces',belements,'Vertices',nodes*1e3,'EdgeColor','yellow','FaceColor','none','LineWidth',2);
    set(gca,'FontSize',12)
    box on
    uananodes = nodes;
    udispnodes = nodes;
    for inode = 1:length(nodes)
        xcoord = nodes(inode,1);
        ycoord = nodes(inode,2);
        [theta,radius] = cart2pol(xcoord,ycoord);
        uradial = (1+nu)*aradius^2*P/E/(bradius^2-aradius^2)*(radius*(1-2*nu)+bradius^2/radius);
        uana = uradial*cos(theta);
        vana = uradial*sin(theta);
        uananodes(inode,1) = uananodes(inode,1) + uana;
        uananodes(inode,2) = uananodes(inode,2) + vana;
    end
    udispnodes(:,1) = udispnodes(:,1) + udisp(1:2:end);
    udispnodes(:,2) = udispnodes(:,2) + udisp(2:2:end);
    patch('Faces',elements,'Vertices',uananodes*1e3,'EdgeColor','black','FaceColor','cyan','FaceAlpha',.5,'LineWidth',1);
    patch('Faces',elements,'Vertices',udispnodes*1e3,'EdgeColor','red','FaceColor','cyan','FaceAlpha',.5,'LineWidth',1);
    patch('Faces',belements,'Vertices',udispnodes*1e3,'EdgeColor','yellow','FaceColor','none','LineWidth',2);
end

%% L2 and H1
L2 = 0;
H1 = 0;
% loop over the elements for error
for iel = 1:numelem
    % get the gauss points and weights
    eType = 4;
    order = 2;
    [w,pt] = gaussptwt (eType,order);
    % get the element connectivity
    econ = elements(iel,:) ;
    % get its coordinates
    nds = nodes(econ,:) ;
    % global index
    gind = [];
    for i = 1:eType
        gind = [gind 2*econ(i)-1 2*econ(i)];
    end
    % elemental solution
    udispele = udisp(gind);

    % loop over the gauss points
    for igp = 1:size(w,2)

        % get the gauss point
        gpt = pt(igp,:);
        % get the shape functions and its derivatives
        [n,dns,dnt] = getShape2d(gpt,eType);
        % to locate the point in the physical space
        xglobal = n*nds(:,1) ;
        yglobal = n*nds(:,2) ;
        % find the jacobian
        jac = [dns;dnt]*nds ;

        % find the derivatives in the physical space
        dnxy = [1 0 0 0; 0 0 0 1; 0 1 1 0]*[inv(jac)*[dns(1) 0 dns(2) 0 dns(3) 0 dns(4) 0; dnt(1) 0 dnt(2) 0 dnt(3) 0 dnt(4) 0]; inv(jac)*[0 dns(1) 0 dns(2) 0 dns(3) 0 dns(4); 0 dnt(1) 0 dnt(2) 0 dnt(3) 0 dnt(4)]] ;
        
        % finite element solution
        uvele = [n(1) 0 n(2) 0 n(3) 0 n(4) 0; 0 n(1) 0 n(2) 0 n(3) 0 n(4)]*udispele;
        % analytical solution
        [theta,radius] = cart2pol(xglobal,yglobal);
        uradial = (1+nu)*aradius^2*P/E/(bradius^2-aradius^2)*(radius*(1-2*nu)+bradius^2/radius);
        uana = uradial*cos(theta);
        vana = uradial*sin(theta);

        % error for u
        uerror = uana - uvele(1);
        verror = vana - uvele(2);
        L2 = L2 + (uerror^2+verror^2)*det(jac)*w(igp);

        % duerror =   - dnxy*udispele;

        H1 = H1 + 0;
    end
end

L2 = sqrt(L2);
H1 = sqrt(H1);

end

%% Gauss Points
function [w,pt] = gaussptwt (eType,order)
if eType == 1 || eType == 2
    if order == 1
        w = 2;
        pt = 0;
    elseif order == 2
        w = [1,1];
        pt = [-1/sqrt(3), 1/sqrt(3)];
    elseif order == 3
        w = [5/9 5/9 8/9];
        pt = [-sqrt(3/5) sqrt(3/5) 0];
    elseif order == 4
        w = [(18-sqrt(30))/36 (18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36];
        pt = [-sqrt(3/7+2/7*sqrt(6/5)) sqrt(3/7+2/7*sqrt(6/5)) -sqrt(3/7-2/7*sqrt(6/5)) sqrt(3/7-2/7*sqrt(6/5))];
    elseif order == 5
        w = [(322-13*sqrt(70))/900 (322-13*sqrt(70))/900 (322+13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225];
        pt = [-1/3*sqrt(5+2*sqrt(10/7)) 1/3*sqrt(5+2*sqrt(10/7)) -1/3*sqrt(5-2*sqrt(10/7)) 1/3*sqrt(5-2*sqrt(10/7)) 0];
    end
elseif eType == 3 || eType == 6
    if order == 1
        w = 0.5;
        pt = [1/3 1/3];
    elseif order == 2
        w = 0;
        pt = [0 0];
    elseif order == 3
        w = 0;
        pt = [0 0];
    end
elseif eType == 4 || eType == 8
    if order == 1
        w = [4];
        pt = [0 0];
    elseif order == 2
        w = [1 1 1 1];
        pt = [1/sqrt(3) 1/sqrt(3); -1/sqrt(3) 1/sqrt(3); -1/sqrt(3) -1/sqrt(3); 1/sqrt(3) -1/sqrt(3)];
    elseif order == 3
        w = [5/9*5/9 5/9*8/9 5/9*5/9 8/9*5/9 8/9*8/9 8/9*5/9 5/9*5/9 5/9*8/9 5/9*5/9];
        pt = [-sqrt(3/5) -sqrt(3/5); -sqrt(3/5) 0; -sqrt(3/5) sqrt(3/5); 0 -sqrt(3/5); 0 0; 0 sqrt(3/5); sqrt(3/5) -sqrt(3/5); sqrt(3/5) 0; sqrt(3/5) sqrt(3/5)];
    end
end
end

%% Get shape functions
function [n,dns,dnt] = getShape2d(pt,eType)
s = pt(1);
t = pt(2);
if eType == 3
    %pt - gauss point
    n = [1-s-t, s, t];
    % derivative is wrt to pt (parameteric coordiante)
    dns = [-1, 1, 0];
    dnt = [-1, 0, 1];
elseif eType == 4
    %pt - gauss point
    n = [0.25*(1-s)*(1-t), 0.25*(1+s)*(1-t), 0.25*(1+s)*(1+t), 0.25*(1-s)*(1+t)];
    % derivative is wrt to pt (parameteric coordiante)
    dns = [-0.25*(1-t), 0.25*(1-t), 0.25*(1+t), -0.25*(1+t)];
    dnt = [-0.25*(1-s), -0.25*(1+s), 0.25*(1+s), 0.25*(1-s)];
elseif eType == 8
    %pt - gauss point
    n = [0.25*(1-s)*(t-1)*(1+s+t), 0.25*(1+s)*(t-1)*(t-s+1), 0.25*(1+s)*(1+t)*(s+t-1), 0.25*(s-1)*(t+1)*(s-t+1), 0.5*(1-s^2)*(1-t), 0.5*(1+s)*(1-t^2), 0.5*(1-s^2)*(1+t), 0.5*(1-s)*(1-t^2)];
    % derivative is wrt to pt (parameteric coordiante)
    dns = [-0.25*(t-1)*(2*s+t), -0.25*(t-1)*(2*s-t), 0.25*(t+1)*(2*s-t), 0.25*(t+1)*(2*s+t), 0.5*(-2*s)*(1-t), 0.5*(1-t^2), 0.5*(-2*s)*(1+t), -0.5*(1-t^2)];
    dnt = [0.25*(1-s)*(2*t+s), 0.25*(1+s)*(2*t-s), 0.25*(1+s)*(2*t+s), -0.25*(s-1)*(2*t-s), -0.5*(1-s^2), -0.5*(1+s)*(2*t), 0.5*(1-s^2), -0.5*(1-s)*(-2*t)];
end
end

function [n,dn] = getShape1d(pt,eType)
if eType == 1
    %pt - gauss point
    n = 1/2*[1 - pt, 1 + pt];
    % derivative is wrt to pt (parameteric coordiante)
    dn = 1/2*[-1, 1];
elseif eType == 2
    %pt - gauss point
    n = 1/2*[-pt+pt^2 2-2*pt^2 pt+pt^2];
    % derivative is wrt to pt (parameteric coordiante)
    dn = 1/2*[-1+2*pt, -4*pt 1+2*pt];
end
end