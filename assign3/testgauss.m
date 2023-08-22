function [] = testgauss(eType,order)
for ieType = eType
    for iorder = order
        [w,pt] = gaussptwt (ieType,iorder);
        disp(w)
        disp(pt)
    end
end
end

function [w,pt] = gaussptwt (eType,order)
order = 2;
if eType == 3 || eType == 6
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