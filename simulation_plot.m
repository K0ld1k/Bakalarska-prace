function simulation_plot(x,y,neighbors,target,Ms,agents,connections)
%  SIMULATION_PLOT Vykreslení simulace
    phi = 0:.1:2*pi;
    for o =1:size(Ms,2)
        k = Ms(3,o) * cos(phi);
        l = Ms(3,o) * sin(phi);
        plot(k + Ms(1,o),l + Ms(2,o),'black')
        fill(k + Ms(1,o),l + Ms(2,o),'black')
        hold on;
    end
    box on
    xlabel("Osa X")
    ylabel("Osa Y")
    plot(x,y,'bluediamond')
    hold on;
    plot(target(1),target(2),'red.','MarkerSize',15)
    hold on;
    if(connections == 1)
    for i = 1:1:agents
        len = length(neighbors{i});
        for j = 2:1:len
            vector = neighbors{i}(j);
            plot([x(i) x(vector)],[y(i) y(vector)],'black')
            hold on;
        end
    end
    end
    hold off;
end