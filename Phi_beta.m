function phi_b = Phi_beta(z,d,h)
%  PHI_BETA Výpočet phi beta (phi_b)
    sig_d = sigma(d);    
    ph_z = (z/sig_d);
    sig_z = (z-sig_d);
    bump_fc = ph(ph_z,h);  
    sig1 = sig_z/sqrt(1 + (sig_z.^2));
    phi_b = bump_fc * (sig1 - 1);       
end

