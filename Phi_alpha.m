function phi_a = Phi_alpha(z,r,d,h)
%  PHI_ALPHA Výpoèet phi alfa (phi_a)
    sig_d = sigma(d);
    sig_r = sigma(r);    
    phi = Phi_z(z - sig_d);
    ph_z = (z/sig_r);
    bump_fc = ph(ph_z,h);    
    phi_a = bump_fc * phi;        
end