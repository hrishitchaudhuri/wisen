clc; clear all; close all;

%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%

% Field Dimensions
x_max = 100;
y_max = 100;

% Sink Coordinates
sink.x = 0.5 * x_max;
sink.y = 0.5 * y_max;

% Number of Nodes
NUM_NODES =100;

% Number of Clusters
k = 7;

% Initial Energy of Each Node
Eo=0.5;
%Radius of mobile sink path
radius = 25;
% Packet Length
packet_length = 500;


% Radio Transmission and Reception Coefficients
ETX=50*10^(-9);
ERX=50*10^(-9);

% Free Space and Multi-Path Model Coefficients
Eelec=50*10^(-9); % units in Joules/bit
Efs= 10*0.0000000001;
Emp= 0.013*0.0000000001;
EDA=5*10^(-9);


% Transmit Amplifier Types %
Eamp=100*10^(-12); % units in Joules/bit/m^2 (amount of energy spent by the amplifier to transmit the bits)


%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%

% Threshold Distance
do = Efs/Emp;

% Creating a random sensor network
figure(1);

for i = 1:1:NUM_NODES
    S(i, 1) = rand(1, 1) * x_max;
    S(i, 2) = rand(1, 1) * y_max;
    
    plot(S(i, 1), S(i, 2), 'red o');
    title 'Wireless Sensor Network';
    xlabel '(m)';
    ylabel '(m)';
    hold on;
end
% cluster = zeros(k, 20000);
% weight = zeros(k, 20000);
% CH_s_in_each_round = zeros(k, 20000);
%%clustering by k-means
[idx, C] = kmeans(S, k);
hold on;

for i = 1:NUM_NODES
    nodes(i).x = S(i, 1);
    nodes(i).y = S(i, 2);
    nodes(i).battery = Eo;
    nodes(i).cluster = (idx(i));
    nodes(i).role=0;
    nodes(i).cond=1;
    nodes(i).rleft=0;
    nodes(i).dst_from_CH=0;
    nodes(i).dist_from_origin=sqrt((nodes(i).x-50)^2+(nodes(i).y-50)^2);
    if (nodes(i).dist_from_origin >= radius)  %%%%Distance of nodes from circumference of the circular path that the mobile sink takes
        nodes(i).dcir = nodes(i).dist_from_origin - radius;
    else
        nodes(i).dcir = radius - nodes(i).dist_from_origin;
    end
end

%%Obtaining array containg node ID's according to cluster numbers
for i = 1:k
    c=1;
    for j = 1:NUM_NODES
        if (nodes(j).cluster == i)
            cluster(c,i)=j;
            c=c+1;            
        end
    end  
    CH_s(i).num_of_nodes = c;
end

operating_nodes = NUM_NODES;
dead_nodes=0;
round = 0;
while dead_nodes < NUM_NODES
    
    %calculating weight
    for i = 1: k
        counter = 1;
        for j = 1 : NUM_NODES
            if(nodes(j).rleft == 0)
                nodes(j).weight = (nodes(j).battery)^2/nodes(j).dcir;
            
            end
            
            if (nodes(j).cluster == i)
                
                weight(counter,i)=nodes(j).weight;
                counter=counter+1;
            end
        end
    end
    
        for j = 1:NUM_NODES
            
                if(nodes(j).rleft > 0)
                    nodes(j).rleft = nodes(j).rleft - 1;
                end
            
        end
    




    %%CH election
    for i = 1:k
        [energy_max,ID]=max(weight(:,i));
        CH_s(i).id=cluster(ID,i);%%Contains ID's of CH's                
        nodes(cluster(ID,i)).role=1;
        nodes(cluster(ID,i)).rleft = 10;
        nodes(cluster(ID,i)).weight = 0;
        CH_s_in_each_round(i,round+1) = CH_s(i).id;
            
        
        
        for j = 1:NUM_NODES
            if (nodes(j).cluster == i)
            nodes(j).dst_from_CH = sqrt((nodes(j).x-nodes(cluster(ID,i)).x)^2 + (nodes(j).y-nodes(cluster(ID,i)).y)^2);%Distance of node from its CH
            end
        end
    end
    %%Energy Consumption during transmission and reception 
    for i = 1:NUM_NODES
        if (nodes(i).cond == 1)
            if (nodes(i).role == 0)
%                 nodes(i).battery = nodes(i).battery - (Eelec*packet_length + Eamp*packet_length*nodes(i).dst_from_CH^2);
              
                    
                if (nodes(i).dst_from_CH < do)
                    nodes(i).battery = nodes(i).battery - (packet_length*Eelec + Efs*packet_length*(nodes(i).dst_from_CH)^2);
                end
                if (nodes(i).dst_from_CH >= do)
                    nodes(i).battery = nodes(i).battery - (packet_length*Eelec + Efs*packet_length*(nodes(i).dst_from_CH)^4);
                end
            end
            if (nodes(i).role == 1)
%                nodes(i).battery = nodes(i).battery - (Eelec+EDA)*packet_length;
                nodes(i).battery = nodes(i).battery - packet_length*CH_s(nodes(i).cluster).num_of_nodes*Eelec;
            end
            
            if (nodes(i).battery <= 0)
                nodes(i).cond = 0;
                dead_nodes = dead_nodes + 1;
                disp("Dead node ID");
                disp(i);
                operating_nodes = operating_nodes - 1;
            end
        
        end
    
    
    end
    %disp(dead_nodes);
     
    round = round + 1;
    op(round)=operating_nodes;
    disp(round);
    disp(dead_nodes);
%     if (dead_nodes == 1 )
%         break;
%     end
    
    for i = 1:NUM_NODES
        if(nodes(i).role == 1)
            nodes(i).role = 0;
        end
    end
     
end

%Plotting Simulation Results "Operating Nodes per Round" %
    figure(2)
    plot(1:round,op(1:round),'-r','Linewidth',2);
    axis([0  20000    0  NUM_NODES]);
    title ({'Updated k_means'; 'Operating Nodes per Round';})
    xlabel 'Rounds ';
    ylabel 'Operational Nodes ';
    hold on;

    



