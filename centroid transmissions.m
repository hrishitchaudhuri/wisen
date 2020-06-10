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
k = 3;

% Initial Energy of Each Node
Eo=0.5;

% Radio Transmission and Reception Coefficients
ETX=50*0.0000001;
ERX=50*0.0000001;

% Free Space and Multi-Path Model Coefficients
Efs=10*0.0000000001;
Emp=0.0013*0.0000000001;

% Number of Rounds
rounds = 100;

%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%

% Threshold Distance
do=sqrt(Efs/Emp);

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
end

wake_clust = 1;
dead_nodes = 0;

for i = 1:1:(k*rounds)
    for i = 1:1:NUM_NODES
        if (nodes(i).state == 1 && nodes(i).cluster == wake_clust)
            
            distance=sqrt((nodes(i).x - nodes(i).sink_x)^2 + (nodes(i).y - nodes(i).sink_y)^2);
            if (distance>do)
                nodes(i).battery = nodes(i).battery - ((ETX)*(4096) + Emp*4096*(distance^4)); 
            end
        
            if (distance<=do)
               nodes(i).battery = nodes(i).battery - ((ETX)*(4096) + Efs*4096*(distance^2)); 
            end
        
            if (nodes(i).battery <= 0)
                dead_nodes = dead_nodes + 1;
                nodes(i).state = 0;
                plot(nodes(i).x, nodes(i).y, 'black x');
            end
            
        end
    end
    
    pause(0.5);
    wake_clust = mod(wake_clust + 1, k);
end