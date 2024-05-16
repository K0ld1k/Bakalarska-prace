function cons_mat = cons_beta_mat(x, y, q_ik, d_b, agents,h)
%  CONS_BETA_MAT Vypocet matice pro consensus agentu beta
cons_mat = zeros(agents,length(q_ik));
for j = 1:1:length(q_ik)
    q_k = q_ik{j};
    for i = 1:1:agents
        z = sqrt((q_k(i,1) - x(i))^2 + (q_k(i,2) - y(i))^2);
        sig = sigma(z);
        sig_db = sigma(d_b);
        cons_mat(i,j) = ph(sig/sig_db,h);
    end
end    
end