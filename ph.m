function bump_fc = ph(z,h)
%  PH Výpoèet bump funkce
% h = 0.2;
if(z >= 0 && z < h)         
    bump_fc = 1;            
elseif(z >= h && z <= 1)
    bump_fc = (1/2) * (1 + cos(pi * ((z - h) / (1 - h))));
else
    bump_fc = 0;
end
end