function [L2,H1] = poisson(nodes,elements,ubcnode,viewplot)
% --- Input Arguments ---
% mesh variables: nodes, elements, ubcnode
% mesh visualization - viewplot
% --- Output Arguments ---
% L2 norm, H1 norm

% mesh details
numelem = size(elements,1); % number of elements
% numbelem = size(belements,1); % number of boundary elements
numnode = size(nodes,1);  % number of nodes
ndof = 1; % number of dof per node

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
    gind = econ;

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
        
        % elemental bilinear and linear form
        kmat(gind,gind) = kmat(gind,gind) + dnxy'*dnxy*w(igp)*det(jac);

        % linear form
        body = 60*yglobal-12*yglobal^2-10-6*xglobal;
        fvec(gind,1) = fvec(gind,1) + n'*body*w(igp)*det(jac);
    end
end

%% Applying Essential Boundary Conditions

% specify the location of essential boundary conditions
bcx = nodes(ubcnode,1);
bcy = nodes(ubcnode,2);
ubcval = bcx.^3+5*bcy.^2-10*bcy.^3+bcy.^4;

TotalDof = 1:ndof*numnode;
activeDof = setdiff(TotalDof,[ubcnode]');

temp = zeros(size(TotalDof,1),1) ;
for in=1:length(ubcnode)
    temp = temp + kmat(:,ubcnode(in))*ubcval(in) ;
end

%modified right hand side
fmod = fvec- temp ;

% modified bilinear form
kmod = kmat(activeDof,activeDof);

% solve
uactive = kmod\fmod(activeDof,1);

% has to be supplemented with the boundary conditions
ufem(ubcnode,1) = ubcval;
ufem(activeDof,1) = uactive ;

%% Plot solution
if viewplot
    meshfig = figure('Name','Mesh','NumberTitle','off');
    grid off;
    axis on;
    xlabel ('x-axis');
    ylabel ('y-axis');
    % xlim([-0.125 1.125]);
    % ylim([-0.125 1.125]);
    daspect([1 1 1]);
    % set(gca, 'color', 'none', 'fontname','calibri');                                              % To remove background
    set(gcf,'units','pixels','position',[500 300 400 400]);                 % To change the size of the figure
    patch('Faces',elements(:,[1 5 2 6 3 7 4 8]),'Vertices',nodes,'FaceVertexCData',ufem,'FaceColor','interp','FaceAlpha',1);
    % patch('Faces',elements(:,[1 5 2 6 3 7 4 8]),'Vertices',nodes,'FaceVertexCData',nodes(:,1).^3+5*nodes(:,2).^2-10*nodes(:,2).^3+nodes(:,2).^4,'FaceColor','interp','FaceAlpha',1);
    hold on
    plot(nodes(:,1),nodes(:,2),'.k');
    plot(nodes(ubcnode,1),nodes(ubcnode,2),'.k','MarkerSize',15);
    set(gca,'FontSize',12)
    colorbar;
    box on
end

%% L2 and H1
L2 = 0;
H1 = 0;
% loop over the elements for error
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
    gind = econ;
    % elemental solution
    ufemelens = ufem(gind);
    
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
        
        % finite element solution
        ufemele = n*ufemelens;
        dufemele = dnxy*ufemelens;
        % analytical solution
        uanaele = xglobal^3+5*yglobal^2-10*yglobal^3+yglobal^4;
        duanaele = [3*xglobal^2; 10*yglobal-30*yglobal^2+4*yglobal^3];

        % error for u
        uerror = (uanaele-ufemele);
        % uerror = (uanaele-ufemele)/uanaele;
        L2 = L2 + uerror'*uerror*det(jac)*w(igp);
        % error for du;
        duerror = (duanaele-dufemele);
        % duerror = (duanaele-dufemele)/norm(duanaele);
        H1 = H1 + duerror'*duerror*det(jac)*w(igp);
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
    dns = [-0.25*(t-1)*(2*s+t), -0.25*(t-1)*(2*s-t), 0.25*(t+1)*(2*s+t), 0.25*(t+1)*(2*s-t), 0.5*(-2*s)*(1-t), 0.5*(1-t^2), 0.5*(-2*s)*(1+t), -0.5*(1-t^2)];
    dnt = [0.25*(1-s)*(2*t+s), 0.25*(1+s)*(2*t-s), 0.25*(1+s)*(2*t+s), -0.25*(s-1)*(2*t-s), -0.5*(1-s^2), -0.5*(1+s)*(2*t), 0.5*(1-s^2), -0.5*(1-s)*(2*t)];
end
end