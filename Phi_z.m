function [phii] = Phi_z(z)
% Vypocet funkce phi
   a = 5;
   b = 5;
   c = abs(a - b)/sqrt(4*a*b);
   z = z + c;
   sig1_z = z/sqrt(1 + (z.^2));        
   phii = (1/2) * (((a + b) * sig1_z) + (a - b)); 
end

