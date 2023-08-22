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
end
end