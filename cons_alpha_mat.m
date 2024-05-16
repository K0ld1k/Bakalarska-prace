function cons_mat = cons_alpha_mat(x, y, neighbors, it_range, agents,h)
%  CONS_ALPHA_MAT Vypocet matice pro consunsus agentu alpha
cons_mat = zeros(agents,agents);
for i = 1:1:agents
    neigh_len = length(neighbors{i});
    for j = 1:1:neigh_len
        n = neighbors{i}(j);
        z = sqrt((x(n) - x(i))^2 + (y(n) - y(i))^2);
        sig = sigma(z);
        sig_it_ran = sigma(it_range);
        cons_mat(i,n) = ph(sig/sig_it_ran,h);
    end
end
end