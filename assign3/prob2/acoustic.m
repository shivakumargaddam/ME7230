function [natf] = acoustic (nodes,elements,c,viewplot)
% --- Input Arguments ---
% mesh variables: nodes, elements, ubcnode
% material constants: c-speed of sound
% mesh visualization: viewplot
% --- Output Arguments ---
% L2 norm, H1 norm



tic
% generate the elements
numelem = length(elements); % number of elements
numnode = length(nodes);  % number of nodes per element
ndof = 1; % number of dof per node

% initialize the global matrices
kmat = zeros(ndof*numnode,ndof*numnode);
mmat = zeros(ndof*numnode,ndof*numnode);
fvec = zeros(ndof*numnode,1);

% loop over the elements
for iel = 1:numelem
    % get the gauss points and weights
    eType = 3;
    order = 1;
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
        kmat(gind,gind) = kmat(gind,gind) + (dnxy'*dnxy)*w(igp)*det(jac) ;
        mmat(gind,gind) = mmat(gind,gind) + n'*(1/c^2)*n*w(igp)*det(jac) ;
    end
end

% solve
natf = sqrt(eig(kmat,mmat));

%% Plot mesh
if viewplot
    meshfig = figure('Name','Mesh','NumberTitle','off');
    grid off;
    axis on;
    patch('Faces',elements,'Vertices',nodes,'FaceColor','cyan','FaceAlpha',1);
    xlabel ('Length [m]');
    ylabel ('Width [m]');
    daspect([1 1 1]);
    set(gca, 'color', 'none', 'fontname','calibri');                       % To remove background and change font
    set(gcf,'units','pixels','position',[500 300 800 600]);                % To change the size of the figure
    box on
end
toc
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
        w = [1/3 1/3 1/3];
        pt = [1/6 2/3; 1/6 1/6; 2/3 1/6];
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