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
    % %pt - gauss point
    % n = 1/2*[-pt+pt^2 2-2*pt^2 pt+pt^2];
    % % derivative is wrt to pt (parameteric coordiante)
    % dn = 1/2*[-1+2*pt, -4*pt 1+2*pt];
end