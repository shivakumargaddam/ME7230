clear all
close all
% clc

%% Shape Functions and Their derivatives w.r.to s,t -- Linear Quadrilateral Element
% Shape Functions
N1  = @(s,t) 0.25*(1-s)*(1-t);
N2  = @(s,t) 0.25*(1+s)*(1-t);
N3  = @(s,t) 0.25*(1+s)*(1+t);
N4  = @(s,t) 0.25*(1-s)*(1+t);
% Derivatives of the shape functions w.r.to 's'
N1s  = @(t) -0.25*(1-t);
N2s  = @(t)  0.25*(1-t);
N3s  = @(t)  0.25*(1+t);
N4s  = @(t) -0.25*(1+t);
% Derivatives of the shape functions w.r.to 't'
N1t  = @(s) -0.25*(1-s);
N2t  = @(s) -0.25*(1+s);
N3t  = @(s)  0.25*(1+s);
N4t  = @(s)  0.25*(1-s);


for order = 1:5
    [gp, gw] = gauss_pw2D(order);                                                  % Gauss points and weights for 2D numerical integration
    Nodes = [0 0; 1 0; 1 1; 0 1];
    NumInt = 0;
    for j = 1:order^2
        % Gauss point and respective weights of each coordinate
        s = gp(j,1);
        t = gp(j,2);
        ws = gw(j,1);
        wt = gw(j,2);
        
        J24 = [N1s(t) N2s(t) N3s(t) N4s(t);
            N1t(s) N2t(s) N3t(s) N4t(s)];
        J = J24*Nodes;
        
        Ni = [N1(s,t) N2(s,t) N3(s,t) N4(s,t)];
        x = Ni*Nodes(:,1);
        y = Ni*Nodes(:,2);
        
%         NumInt = NumInt + (ws*wt)*(1+x^2+y^2-2*y^3)*det(J); % Problem 1
%         NumInt = NumInt + (ws*wt)*((1-x)*(1-y)*x*(1-y))*det(J); % Problem 4
%         NumInt = NumInt + (ws*wt)*((1-2*x)*(1+y^2-2*y))*det(J); % Problem 5a
        NumInt = NumInt + (ws*wt)*((x-x^2)*(2*y-2))*det(J); % Problem 5b
    end
    
    sprintf('no of points = %d, value = %0.15f',order^2,NumInt)
end

function [gp, gw] = gauss_pw2D (order)
%% HELP - gauss_pw2D
%   -- generates gauss points and weights for 2D numerical integration
%   -- input 'N' gives, 'N^2' no of points and weights
%   -- N: even numbers only [No of gauss points in each dimension]

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

x = pt';
wx = w';
y = pt';
wy = w';
gp = zeros(order^2, 2);
gw = zeros(order^2, 2);

row = 1;
for i = x'
    for j = y'
        gp(row,:) = [i j];
        row = row+1;
    end
end

row = 1;
for i = wx'
    for j = wy'
        gw(row,:) = [i j];
        row = row+1;
    end
end

end