clear all
close all
clc

% nodal coordiantes
node = [0,0; 1,1; 0,1; 1,0; 0.5,0.5] ;
element = [1 4 5; 3 5 2; 1 5 3; 5 4 2];


% node = [0,0; 1,1; 0,1; 1,0; 0.5,0.5;0 0.5;1 0.5] ;
% element = [1 4 5;3 5 2;6 5 3;1 5 6;5 4 7;5 7 2];

eType = 3;
k = 1;

% generate the elements
numelem = length(element); % number of elements
numnode = size(node,1);
nnode = length(node);  % number of nodes per element
ndof = 1; % number of dof per node


% initialize the global matrices
kmat = zeros(ndof*nnode,ndof*nnode);
fvec = zeros(ndof*nnode,1);

% gauss points
order = 1;
if order == 1
    w = 0.5;
    pt = [1/3 1/3];
elseif order == 2

elseif order == 3

end

% loop over the elements
for iel = 1:numelem
    
    % get the element connectivity
    econ = element(iel,:) ;
    
    % get its coordinates
    nds = node(econ,:) ;
    
    % global index
    gind = econ ;
    
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
        dnxy = inv(jac)*[dns;dnt] ;
        
        % elemental bilinear and linear form
        kmat(gind,gind) = kmat(gind,gind) + dnxy'*k*dnxy*w(igp)*det(jac);
        
        % linear form
        b = 0;
        fvec(gind,1) = fvec(gind,1) + n'*b*w(igp)*det(jac);
        
    end
end

% specify the location of boundary conditions
bcdof = [1 2 3 4];
bcval = [10 2 10 2];
% 
% bcdof = [1 2 3 4 6 7];
% bcval = [10 2 10 2 10 2];

TotalDof = 1:ndof*numnode ;
activeDof = setdiff(TotalDof,bcdof);


temp = zeros(size(TotalDof,1),1) ;
for in=1:length(bcdof)
    temp = temp + kmat(:,bcdof(in))*bcval(in) ;
end

%modified right hand side
fmod = fvec- temp ;

% modified bilinear form
kmod = kmat(activeDof,activeDof) ;

% solve
u = kmod\fmod(activeDof,1)

% has to be supplemented with the boundary conditions
udisp = zeros(ndof*numnode,1) ;
udisp(bcdof) = bcval ;
udisp(activeDof) = u ;

udisp

% figure; hold on;
% plot(node,udisp,'rs-')
% 
% % % plot the analytical solution
% % x = linspace(xo,x1,100);
% % for i=1:length(x)
% %     uana(i,1) = 2/3*x(i)^3 + 34/3*x(i);
% % end
% % 
% % plot(x,uana,'b-')
% 
% % loop over the elements for error
% L2 = 0;
% H1 = 0;
% for iel = 1:numelem
% 
%     % get the element connectivity
%     econ = element(iel,:) ;
% 
%     % get its coordinates
%     nds = node(econ) ;
% 
%     % global index
%     gind = econ ;
%     % elemental solution
%     uele = udisp(gind);
% 
%     % l = 0 ;
%     ele = zeros(2,2);
%     % loop over the gauss points
%     for igp = 1:size(w,2)
% 
%         % get the gauss point
%         gpt = pt(igp);
% 
%         % get the shape functions and its derivatives
%         [n,dn] = getShape1d(gpt);
% 
%         % to locate the point in the physical space
%         xglobal = n*nds' ;
% 
%         % plot the gauss point in the physical space
%         %plot(xglobal(1),0,'bs');
% 
%         % find the jacobian
%         jac = dn*nds' ;
% 
%         % find the derivatives in the physical space
%         dnx = inv(jac)*dn ;
% 
%         % error for u
%         uerror = 2/3*xglobal^3 + 34/3*xglobal - n*uele;
%         L2 = L2 + uerror*uerror*det(jac)*w(igp);
%         duerror = 2*xglobal^2 + 34/3 - dnx*uele;
%         H1 = H1 + duerror*duerror*det(jac)*w(igp);
% 
%         % error for u'
% 
%     end
% end
% 
% L2 = sqrt(L2)
% H1 = sqrt(H1)
% h = node(2) - node(1)
