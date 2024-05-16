function beta_plot(q_ik,r_b,x,y,agents)
    % Vykreslení vazeb mezi alpha a beta agenty
     for i = 1:1:length(q_ik)
        q_k = q_ik{i};
        for j = 1:1:agents
            distance = sqrt((q_k(j,1) - x(j))^2 + (q_k(j,2) - y(j))^2);
            if(distance <= r_b)
                plot([x(j) q_k(j,1)],[y(j) q_k(j,2)],'red')
                hold on;
                plot(q_k(j,1),q_k(j,2),'greendiamond')
            end
        end
    end
end

