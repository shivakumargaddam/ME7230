clear all
close all
% clc

%% Shape Functions and Their derivatives
% Shape Functions
N1  = @(s) 0.5*(1-s);
N2  = @(s) 0.5*(1+s);
% Derivatives of the shape functions w.r.to 's'
N1s  = -0.5;
N2s  =  0.5;



%%
for order = 1:5
    %% Gauss points 2D
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
    %% Gauss Integration
    Nodes = [0 1];
    NumInt = 0;
    for j = 1:order
        % Gauss point and respective weights of each coordinate
        s = pt(j);
        ws = w(j);

        J24 = [N1s N2s];
        J = J24*Nodes';
        
        Ni = [N1(s) N2(s)];
        x = Ni*Nodes';
        
        NumInt = NumInt + (ws)*(1/x)*det(J); % Problem 3
    end
    
    sprintf('no of points = %d, value = %0.15f',order,NumInt)
end
