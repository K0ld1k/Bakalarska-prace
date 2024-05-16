function [neighbors] = gen_neighbors(x, y,neighbors, agents, it_range)
% GEN_NEIGHBORS Funkce d�v� po�et soused� ka�d�ho agenta
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
        % Z�sk�n� index� soused� v�etn� agenta samotn�ho
        all_neighbors = find(neigh(i,:) < it_range);
        
        % Ujist�me se, �e agent s�m se sebou je na prvn�m m�st�
        all_neighbors(all_neighbors == i) = [];
        neighbors{i} = [i, all_neighbors];
    end
end