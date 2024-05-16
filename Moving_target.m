clear all
close all
clc
%% INICIALIZACE
% Pocet agentu, dimenze, vzdalenost, skalovani, reakcni vzdalenost
agents = 30; %<n>
dim = 2; %<dim>
d = 17; %<d>
K = 1.2; %<k>
r = K*d; %<r>

% Rozmisteni agentu a jejich sousedi
a = 50;
b = 150;
aa = -50;
bb = 50;
x = a + (b-a).*rand(agents,1)
y = aa + (bb-aa).*rand(agents,1)
neighbors = {};
[neighbors] = gen_neighbors(x,y,neighbors,agents,r);

% Prekazka a cil
Ms = [200 200 200 320 320 320 450 450 450;   %x 
      -75   0  75 -75   0  75 -75   0  75;   %y
       20  20  20  20  20  20  20  20  20];  %radius 

target = [550,0];
start = 1;
connections = 0;
tolerance = 0.1;

%Promenne pro vypocet vzdalenosti od prekazky
d_b = 0.7 * d;
r_b = 1.2 * d_b;

% Vykres pocatecni polohy
figure
simulation_plot(x,y,neighbors,target,Ms,agents,connections)
axis([0 650 -200 200])

% Parametry rizeni
epsilon = 0.1; %0<Epsilon<1
h_alpha = 0.2;
h_beta = 0.7;
%c>0
c1_alpha = 30;
c2_alpha = 2 * sqrt(c1_alpha);

c1beta = 600;
c2beta = 2 * sqrt(c1beta);

c1_gama = 20;
c2_gama = 2 * sqrt(c1_gama);
% Doba simulace
%----------------------
d_t = 0.03; %(0,01 - 0,03)
t = 0:d_t:60;
%----------------------
%Posun cile
posun_plus = d_t*100;
posun_minus = d_t*100*(-1);
%% ALGORITMUS
% Priprava pohybu agentu
old_x = x;
old_y = y;
% Rychlost agentu
nodes_vel1 = zeros(agents,dim);
% Rizeni
u_i = zeros(agents,dim);
% Trajektorie
x_iter = zeros(agents,length(t));
y_iter = zeros(agents,length(t));
% Rizeni-zasahy
ua_iter = zeros(2,length(t));
ub_iter = zeros(2,length(t));
uc_iter = zeros(2,length(t));
% Trajektorie stredu shluku
center_iteration = zeros(length(t),2);
% Konektivita
con_iter = zeros(length(t),1);

% Jednotlive rychlosti agentu
velocity_it = zeros(length(t),agents);
% Iterace pro graf rychlosti
iter_number = ones(length(t),1);

for iteration = 1:length(t)
    % Vypocet a ukladani rychlosti
    nodes_vel1(:,1) = (x - old_x)/d_t;
    nodes_vel1(:,2) = (y - old_y)/d_t;

    % Kvadrat rychlosti pro graf + iterace
    if(iteration > 1)
       velocity_it(iteration,:) = sqrt((nodes_vel1(:,1).^2)+(nodes_vel1(:,2).^2));
       iter_number(iteration) = iteration;
    end

    % Ukladani predchozi polohy agentu
    old_x = x;
    old_y = y;
    
    % Vypocet stredu shluku 
    xy_center = zeros(1,2);
    xy_center(1) = sum(x)/agents;
    xy_center(2) = sum(y)/agents;

    % Ulozeni pro ukazku prubehu trajektorie
    center_iteration(iteration,:) = xy_center;

    %Vypocet novych sousedu
    [neighbors] = gen_neighbors(x,y,neighbors,agents,r);
    
    %Vypocet matice sousednosti + kontrola konektivity
    %<neighbors, agents>
    adj_mat = zeros(agents,agents);
    for i = 1:1:agents
        neigh_len1 = length(neighbors{i});
        for j = 2:1:neigh_len1
            n = neighbors{i}(j);
            adj_mat(i,n) = 1;
        end
    end

    % Ulozeni konektivity jednotlivych iteraci
    con_iter(iteration) = rank(adj_mat);

    %-------------------------------------------------------------------------------------------------------
    % Vypocet pozice a rychlosti beta agenta
    %<x, y, obstacle_tran, radius, agents, agent_vel1>
    q_k = zeros(agents,2);
    q_ik = {size(Ms,2)};
    p_ik = {size(Ms,2)};
    I = eye(agents);
    eucld_dist_obst = zeros(agents,1);
    for j = 1:1:size(Ms,2)
        obstacle = Ms(:,j); 
        obs_tran = obstacle(1:2);
        obs_radius = obstacle(3);

        for i = 1:1:agents
            eucld_dist_obst(i) = sqrt((obs_tran(1) - x(i))^2 + (obs_tran(2) - y(i))^2);
        end    
        mu = obs_radius./eucld_dist_obst;
        a_k = [(x - obs_tran(1))./eucld_dist_obst,(y - obs_tran(2))./eucld_dist_obst];
        P = I - (a_k * a_k');
        p_k = P * nodes_vel1;
        for i = 1:1:agents
            p_k = mu(i) * p_k;
            q_k(i,:) = (mu(i) * [x(i),y(i)]) + ((1 - mu(i)) * obs_tran');
        end
        q_ik{j} = q_k;
        p_ik{j} = p_k;
    end
    %-------------------------------------------------------------------------------------------------------
    % Gradient base term alpha
    %<x, y, distance, it_range, neighbors, agents>
    grad_alpha = zeros(agents,2);
    for i = 1:1:agents
        neigh_len2 = length(neighbors{i});
            nij_val = zeros(agents,2);
            for j = 2:1:neigh_len2
                n1 = neighbors{i}(j);
                eucld_dist_a = sqrt((x(n1) - x(i))^2 + (y(n1) - y(i))^2);
                nij_val(n1,1) = (x(n1) - x(i))/sqrt(1 + (epsilon * (eucld_dist_a^2)));
                nij_val(n1,2) = (y(n1) - y(i))/sqrt(1 + (epsilon * (eucld_dist_a^2)));
            end
        for k = 2:1:neigh_len2
            n2 = neighbors{i}(k);
            z = sqrt((x(n2) - x(i))^2 + (y(n2) - y(i))^2);
            sigma1 = sigma(z);                
            grad_alpha(i,:) = grad_alpha(i,:) + (Phi_alpha(sigma1,r,d,h_alpha) * nij_val(n2,:)); 
        end
    end
    
    % Consensus base term alpha
    %<nodes_vel1, neighbours, con_mat_a, agents>
    cons_aplha = zeros(agents,2);
    con_mat_a = cons_alpha_mat(x,y,neighbors,r,agents,h_alpha);
    for i = 1:1:agents
        neigh_len3 = length(neighbors{i});
        for j = 2:1:neigh_len3
            n3 = neighbors{i}(j);
            cons_aplha(i,:) = cons_aplha(i,:) + ((con_mat_a(i,n3) * (nodes_vel1(j,:) - nodes_vel1(i,:)))); 
        end
    end
    %-------------------------------------------------------------------------------------------------------
    % Gradient base term (beta agent)
    %<x, y, q_k, r_beta, d_beta, agents>
    grad_beta1 = zeros(agents,2);
    grad_beta = zeros(agents,2);
    nik_val = zeros(agents,2);
    for k = 1:1:length(q_ik)
        q_1 = q_ik{k};
            for i = 1:1:agents
                eucld_dist_b = sqrt((q_1(i,1) - x(i))^2 + (q_1(i,2) - y(i))^2);
                nik_val(i,1) = (q_1(i,1) - x(i))/sqrt(1 + (epsilon * (eucld_dist_b^2)));
                nik_val(i,2) = (q_1(i,2) - y(i))/sqrt(1 + (epsilon * (eucld_dist_b^2)));
            end
        for j = 1:1:agents
            z = sqrt((x(j) - q_1(j,1))^2 + (y(j) - q_1(j,2))^2);
            sig1 = sigma(z);
            grad_beta1(j,:) = Phi_beta(sig1,d_b,h_beta) * nik_val(j,:);
        end
        grad_beta = grad_beta + grad_beta1;
    end
    % Consensus base term (beta agent)
    %<node_vel1, p_k, con_mat_b, agents>
    cons_beta = zeros(agents,2);
    con_mat_b = cons_beta_mat(x, y, q_ik, d_b, agents,h_beta);
    for j = 1:1:length(p_ik)
        p_k = p_ik{j};
        for i = 1:1:agents
            cons_beta(i,:) = cons_beta(i,:) + con_mat_b(i,j) * (p_k(i,:) - nodes_vel1(i,:));
        end
    end
    %-------------------------------------------------------------------------------------------------------
    % Upraveni gradientu a consensu o doporucene parametry
    grad_alpha = c1_alpha * grad_alpha;
    cons_aplha = c2_alpha * cons_aplha;
    grad_beta = c1beta * grad_beta;
    cons_beta = c2beta * cons_beta;
    
    % Vypocet zmeny pozice vuci cilu
    position_dif = zeros(agents,2);
    position_dif(:,1) = (x - target(1));
    position_dif(:,2) = (y - target(2));
    
    % Vypocet rizeni
    u_aplha = grad_alpha + cons_aplha;
    u_beta = grad_beta + cons_beta;
    sig_gamma = sigma_1(position_dif);
    u_gama1 = (c1_gama * sig_gamma);
    u_gama2 = (c2_gama * nodes_vel1);
    u_i = u_aplha + u_beta - u_gama1 - u_gama2;
    
    % Zmena pozic agentu
    x = old_x + (d_t * nodes_vel1(:,1)) + (((d_t^2)/2) * u_i(:,1));
    y = old_y + (d_t * nodes_vel1(:,2)) + (((d_t^2)/2) * u_i(:,2));
    
    % Ukladani hodnot kazde iterace
    x_iter(:,iteration) = x;
    y_iter(:,iteration) = y;
    ua_iter(:,iteration) = u_aplha(5,:);
    ub_iter(:,iteration) = u_beta(5,:);
    ug = - u_gama1 - u_gama2;
    uc_iter(:,iteration) = ug(1,:);

    % Posun po X
    % target(1) = target(1)+d_t*10;

    % Posun po Y
    if(start ~= 1)
       target(2) = target(2)+posun; 
    end
    if(target(2) == 0 && start == 1)
        target(2) = target(2)+posun_plus;
        posun = posun_plus;
        start = 0;
    elseif(target(2) >= 150)
       posun = posun_minus;
    elseif(target(2) <= -150)
       posun = posun_plus;
    end
    if(rank(adj_mat) == agents && iteration>500)
        posun_plus = 0;
        posun_minus = 0;
        posun = 0;
    end

    hold off;
    simulation_plot(x,y,neighbors,target,Ms,agents,connections)
    hold on
    if(connections == 1)
    beta_plot(q_ik,r_b,x,y,agents)
    end
    axis([0 650 -200 200])
    hold off;
    drawnow;

    if(xy_center(1) > target(1)-tolerance && xy_center(1) < target(1)+tolerance && xy_center(2) > target(2)-tolerance && xy_center(2) < target(2)+tolerance && rank(adj_mat) == agents)
        sim_end = iteration;
        break
    else
        sim_end = iteration;
    end

end

%% VYSLEDKY
pocet_iteraci_simulace = length(t)
konec_simulace = sim_end
delka_simulace = length(t)-(length(t)-sim_end)
x_iter = x_iter';
y_iter = y_iter';
%% Prubeh bump funkce
% figure
% values = 0:0.01:2;
% ph_values = zeros(1,length(values));
% for i = 1:length(values)
%     ph_values(i) = ph(values(i),h_alpha); 
% end
% plot(values,ph_values)
% hold on
% for i = 1:length(values)
%     ph_values(i) = ph(values(i),h_beta); 
% end
% plot(values,ph_values)
% grid on
% box on
% title("Nárazová funkce")
% xlabel("z")
% ylabel("\rho_h(z)")
% legend("Nárazová funkce (h=0.2)","Nárazová funkce (h=0.7)")
% 
% %% Graf akcni funkce 
% figure
% values1 = 1:0.1:40;
% phi_alpha_values = zeros(1,length(values1));
% for i = 1:length(values1)
%     phi_alpha_values(i) = Phi_alpha(values1(i),1.2*10,10,h_alpha); 
% end
% plot(values1,phi_alpha_values)
% hold on
% for i = 1:length(values1)
%     phi_alpha_values(i) = Phi_alpha(values1(i),1.2*10,10,h_beta); 
% end
% plot(values1,phi_alpha_values)
% grid on
% box on
% title("Akční funkce")
% xlabel("z")
% ylabel("\phi_\alpha(z)")
% legend("Akční funkce (h=0.2)","Akční funkce (h=0.7)")

%% Graf pro rychlost
figure
for i = 1:1:agents
    plot(iter_number(1:sim_end,1), velocity_it(1:sim_end,i));    
    title('Rychlosti jednotlivých agentů');
    hold on;
    grid on;
    box on;
    set(gcf,'color','w');
end
xlabel("Iterace")
ylabel("Rychlost")

%% Graf konektivity
figure
plot(iter_number(1:sim_end,1), con_iter(1:sim_end,1));
grid on
box on
set(gcf,'color','w');
title('Průběh konektivity shluku');
xlabel("Iterace")
ylabel("Stupeň konektivity")

%% Graf trajektorie
figure
scatter(x_iter(sim_end,:),y_iter(sim_end,:),'bluediamond');
hold on
for i = 1:1:agents
    plot(x_iter(1:sim_end,i),y_iter(1:sim_end,i),'b');
    hold on
     phi = 0:.1:2*pi;
    for o =1:size(Ms,2)
        k = Ms(3,o) * cos(phi);
        l = Ms(3,o) * sin(phi);
        plot(k + Ms(1,o),l + Ms(2,o),'black')
        fill(k + Ms(1,o),l + Ms(2,o),'black')
        hold on;
    end
    box on
    grid on
    set(gcf,'color','w');
end
title('Trajektorie jednotlivých \alpha-agentů')
xlabel("Osa X")
ylabel("Osa Y")

%% Graf stredu shluku
figure
plot(center_iteration(1:sim_end,1),center_iteration(1:sim_end,2),'b-');
hold on
plot(target(1),target(2),'k.');
title('Trajektorie středu shluku')
xlabel("Osa X")
ylabel("Osa Y")
hold on
grid on
box on
set(gcf,'color','w');

%% Graf rizeni 
figure
hold on
plot(iter_number(1:sim_end),ua_iter(1,1:sim_end))
plot(iter_number(1:sim_end),ua_iter(2,1:sim_end))

plot(iter_number(1:sim_end),ub_iter(1,1:sim_end))
plot(iter_number(1:sim_end),ub_iter(2,1:sim_end))

plot(iter_number(1:sim_end),uc_iter(1,1:sim_end))
plot(iter_number(1:sim_end),uc_iter(2,1:sim_end))
title("Řízení agentů")
xlabel("Iterace")
ylabel("Vstup řízení")
legend("u_\alpha(x)","u_\alpha(y)", "u_\beta(x)","u_\beta(y)", "u_\gamma(x)", "u_\gamma(y)")
grid on
box on
set(gcf,'color','w');

%% Tvar vysledneho shluku
figure
plot(digraph(adj_mat))
title("Digraf")
set(gcf,'color','w');

    
    
    
    
    
    
    
