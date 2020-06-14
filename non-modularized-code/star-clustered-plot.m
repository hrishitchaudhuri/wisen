clc; clear; close;

%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%

% Field Dimensions
x_max = 100;
y_max = 100;

% Sink Coordinates
sink.x = 0.5 * x_max;
sink.y = 0.5 * y_max;

% Number of Nodes
NUM_NODES = 100;

% Number of Clusters
k = 5;

% Initial Energy of Each Node
Eo=0.5;

% Packet Length
packet_length = 500;

% Radio Transmission and Reception Coefficients
ETX=50*0.0000001;
ERX=50*0.0000001;

% Free Space and Multi-Path Model Coefficients
Efs= 10*0.0000000001;
Emp= 0.013*0.0000000001;

% Number of Rounds
rounds = 100;

%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%

% Threshold Distance
do = Efs/Emp;

% Creating a random sensor network
figure(1);

for i = 1:1:NUM_NODES
    S(i, 1) = rand(1, 1) * x_max;
    S(i, 2) = rand(1, 1) * y_max;
    
    plot(S(i, 1), S(i, 2), 'red o');
    hold on;
end

% disp(S);

% plot(sink.x, sink.y, 'black o');
% hold on;

[idx, C] = kmeans(S, k);

% disp(idx);
% disp(C);

plot(C(:,1), C(:,2), 'k*', 'MarkerSize', 15);
hold on;

% nodes = zeros;
for i = 1:1:NUM_NODES
    nodes(i).x = S(i, 1);
    nodes(i).y = S(i, 2);
    nodes(i).battery = Eo;
    nodes(i).cluster = (idx(i) - 1);
    nodes(i).sink_x = C(idx(i), 1);
    nodes(i).sink_y = C(idx(i), 2);
    nodes(i).state = 1;
    nodes(i).distance = (nodes(i).x - nodes(i).sink_x)^2 + (nodes(i).y - nodes(i).sink_y)^2;
end

% distances = zeros;
% for i = 0:1:k-1
%     for j = 1:1:NUM_NODES
%         if nodes(j).cluster == i && nodes(i).distance < do
%         end
%     end
% end

wake_clust = 1;
dead_nodes = 0;
plot_idx = 1;
round_stat = zeros;
dead_stats = zeros;

while dead_nodes < NUM_NODES
    for j = 1:1:NUM_NODES
        if (nodes(j).state == 1 && nodes(j).cluster == wake_clust)
            
            if (nodes(j).distance > do)
                nodes(j).battery = nodes(j).battery - ((ETX)*(packet_length) + Emp*packet_length*(nodes(j).distance^2)); 
            end
        
            if (nodes(j).distance <= do)
               nodes(j).battery = nodes(j).battery - ((ETX)*(packet_length) + Efs*packet_length*(nodes(j).distance)); 
            end
        
            if (nodes(j).battery <= 0)
                dead_nodes = dead_nodes + 1;
                nodes(j).state = 0;
                
                figure(1);
                plot(nodes(j).x, nodes(j).y, 'black x');
                hold on;
            end
            
        end
    end
    
    pause(0.05);
    wake_clust = mod(wake_clust + 1, k);
    
    if (mod(i, k) == 0)
        round_stat(plot_idx) = plot_idx;
        dead_stats(plot_idx) = dead_nodes;
        plot_idx = plot_idx + 1;
    end
    
    % if (dead_nodes == NUM_NODES)
    %     break;
    % end
end

figure(2);
plot(round_stat, dead_stats);
hold on;
