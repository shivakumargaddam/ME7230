function [n,dn] = getShape1d(pt,oelement)

if oelement == 1
    %pt - gauss point
    n = 1/2*[1 - pt, 1 + pt];
    % derivative is wrt to pt (parameteric coordiante)
    dn = 1/2*[-1, 1];

elseif oelement == 2
    %pt - gauss point
    n = 1/2*[-pt+pt^2 2-2*pt^2 pt+pt^2];
    % derivative is wrt to pt (parameteric coordiante)
    dn = 1/2*[-1+2*pt, -4*pt 1+2*pt];
end