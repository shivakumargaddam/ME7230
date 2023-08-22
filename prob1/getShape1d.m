function [n,dn] = getShape1d(pt)

%pt - gauss point

n = 1/2*[1 - pt, 1 + pt];

% derivative is wrt to pt (parameteric coordiante
dn = 1/2*[-1, 1];