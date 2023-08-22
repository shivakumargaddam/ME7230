clear all
close all
% clc

%% Shape Functions and Their derivatives
% Shape Functions
N1  = @(s,t) (1-s-t);
N2  = @(s,t) s;
N3  = @(s,t) t;
% Derivatives of the shape functions w.r.to 's'
N1s  = @(t) -1;
N2s  = @(t)  1;
N3s  = @(t)  0;
% Derivatives of the shape functions w.r.to 't'
N1t  = @(s) -1;
N2t  = @(s)  0;
N3t  = @(s)  1;



for order = 1:2
    %% Gauss points 2D
    if order == 1
        gw = 1/2;
        gp = [1/3 1/3];
    elseif order == 2
        gw = [1/6 1/6 1/6]';
        gp = [1/2 1/2; 1/2 0; 0 1/2];
    end
    %% Numerical Integration
    Nodes = [2 -3; 5 -3; 2 0]; % problem 2
    NumInt = 0;
    for j = 1:size(gp,1)
        % Gauss point and respective weights of each coordinate
        s = gp(j,1);
        t = gp(j,2);
        w = gw(j,1);
        
        J24 = [N1s(t) N2s(t) N3s(t);
               N1t(s) N2t(s) N3t(s)];
        J = J24*Nodes;
        
        Ni = [N1(s,t) N2(s,t) N3(s,t)];
        x = Ni*Nodes(:,1);
        y = Ni*Nodes(:,2);
        
        NumInt = NumInt + (w)*(x*y)*det(J); % Problem 2
    end
    
    sprintf('no of point = %d, value = %0.15f',size(gp,1),NumInt)
end

