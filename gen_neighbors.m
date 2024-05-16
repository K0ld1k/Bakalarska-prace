function [neighbors] = gen_neighbors(x, y,neighbors, agents, it_range)
% GEN_NEIGHBORS Funkce dává poèet sousedù každého agenta
    neigh = zeros(agents, agents);
    for i = 1:agents
        for j = 1:agents
            if i == j
                neigh(i,j) = 0; 
            else
                eucld_dist = sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2);
                neigh(i,j) = eucld_dist;
            end
        end
    end
    
    for i = 1:agents
        % Získání indexù sousedù vèetnì agenta samotného
        all_neighbors = find(neigh(i,:) < it_range);
        
        % Ujistíme se, že agent sám se sebou je na prvním místì
        all_neighbors(all_neighbors == i) = [];
        neighbors{i} = [i, all_neighbors];
    end
end