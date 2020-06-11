clc;clear;close;

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
global k;
k = 3;

% Initial Energy of Each Node
Eo=0.5;

% Radio Transmission and Reception Coefficients
global ETX;
global ERX;
ETX=50*0.0000001;
ERX=50*0.0000001;

% Free Space and Multi-Path Model Coefficients
global Efs;
global Emp;
Efs=10*0.0000000001;
Emp=0.0013*0.0000000001;

% Number of Rounds
global rounds;
rounds = 100;

%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%

% Threshold Distance
global do;
do=sqrt(Efs/Emp);

% Creating a random sensor network
figure(1);

% for i = 1:1:NUM_NODES
%     S(i, 1) = rand(1, 1) * x_max;
%     S(i, 2) = rand(1, 1) * y_max;
%     
%     plot(S(i, 1), S(i, 2), 'red o');
%     hold on;
% end
[nodes,pos_nodes] = initialise(NUM_NODES);
% disp(S);

% plot(sink.x, sink.y, 'black o');
% hold on;

[idx, C] = kmeans(pos_nodes, k);

% disp(idx);
% disp(C);

plot(C(:,1), C(:,2), 'k*', 'MarkerSize', 15);
hold on;

% nodes = zeros;
% for i = 1:1:NUM_NODES
%     nodes(i).x = S(i, 1);
%     nodes(i).y = S(i, 2);
%     nodes(i).battery = Eo;
%     nodes(i).cluster = (idx(i) - 1);
%     nodes(i).sink_x = C(idx(i), 1);
%     nodes(i).sink_y = C(idx(i), 2);
%     nodes(i).state = 1;
% end
[nodes] = assign_CHs(nodes,NUM_NODES,C,idx);

wake_clust = 1;
dead_nodes = 0;
plot_idx = 1;
round_stat = zeros;
dead_stats = zeros;

%  for i = 1:1:(k*rounds)
%     for j = 1:1:NUM_NODES
%         if (nodes(j).state == 1 && nodes(j).cluster == wake_clust)
%             
%             distance=sqrt((nodes(j).x - nodes(j).sink_x)^2 + (nodes(j).y - nodes(j).sink_y)^2);
%             if (distance>do)
%                 nodes(j).battery = nodes(j).battery - ((ETX)*(4096) + Emp*4096*(distance^4)); 
%             end
%         
%             if (distance<=do)
%                nodes(j).battery = nodes(j).battery - ((ETX)*(4096) + Efs*4096*(distance^2)); 
%             end
%         
%             if (nodes(j).battery <= 0)
%                 dead_nodes = dead_nodes + 1;
%                 nodes(j).state = 0;
%                 plot(nodes(j).x, nodes(j).y, 'black x');
%             end
%             
%         end
%      end
% 
%     
% %     pause(0.5);
%     wake_clust = mod(wake_clust + 1, k);
%     
%     if (mod(i, k) == 0)
%         round_stat(plot_idx) = plot_idx;
%         dead_stats(plot_idx) = dead_nodes;
%         plot_idx = plot_idx + 1;
%     end
%     
%     if (dead_nodes == NUM_NODES)
%         break;
%     end
% end

[nodes,dead_nodes,round_stat,dead_stats,wake_clust] = energy_model(nodes,NUM_NODES,wake_clust,plot_idx,dead_nodes);

figure(2);
plot(round_stat, dead_stats);
hold on;



%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%
function[nodes,pos_nodes] = initialise(no_of_nodes)
    hold on;
     for i=1:1:no_of_nodes
        nodes(i) = Sensor_Node;
        nodes(i)= Add_Positions(nodes(i));
        pos_nodes(i,1)=nodes(i).x;
        pos_nodes(i,2)=nodes(i).y;
        plot(nodes(i).x,nodes(i).y,'red o');
     end
     nodes(no_of_nodes+1).x = 50;
     nodes(no_of_nodes+1).y = 50;
%      plot(S(no_of_nodes+1).xd,S(no_of_nodes+1).yd,'black x');
end

function[nodes] = assign_CHs(nodes,NUM_NODES,C,idx)
    for i = 1:1:NUM_NODES
    nodes(i).cluster = (idx(i) - 1);
    nodes(i).sink_x = C(idx(i), 1);
    nodes(i).sink_y = C(idx(i), 2);
    end
    
end


function[nodes,dead_nodes,round_stat,dead_stats,wake_clust] = energy_model(nodes,NUM_NODES,wake_clust,plot_idx,dead_nodes)
    hold on;
    global k;
    global rounds;
    global do;
    global ETX;
    global ERX;
    global Efs;
    global Emp;
    
    for i = 1:1:(k*rounds)
    for j = 1:1:NUM_NODES
        if (nodes(j).state == 1 && nodes(j).cluster == wake_clust)
            
            distance=sqrt((nodes(j).x - nodes(j).sink_x)^2 + (nodes(j).y - nodes(j).sink_y)^2);
            if (distance>do)
                nodes(j).battery = nodes(j).battery - ((ETX)*(4096) + Emp*4096*(distance^4)); 
            end
        
            if (distance<=do)
               nodes(j).battery = nodes(j).battery - ((ETX)*(4096) + Efs*4096*(distance^2)); 
            end
        
            if (nodes(j).battery <= 0)
                dead_nodes = dead_nodes + 1;
                nodes(j).state = 0;
                plot(nodes(j).x, nodes(j).y, 'black x');
            end
            
        end
    end
    
     pause(0.5);
    wake_clust = mod(wake_clust + 1, k);
    
    if (mod(i, k) == 0)
        round_stat(plot_idx) = plot_idx;
        dead_stats(plot_idx) = dead_nodes;
        plot_idx = plot_idx + 1;
    end
    
    if (dead_nodes == NUM_NODES)
        break;
    end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%