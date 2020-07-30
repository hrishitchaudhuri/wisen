clc; clear all; close all;

%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%

% Field Dimensions
x_max = 100;
y_max = 100;

% Sink Coordinates
sinkx = 0.5 * x_max;
sinky = 0.5 * y_max;

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

threshold = 70;

% Free Space and Multi-Path Model Coefficients
Eelec=50*10^(-9); % units in Joules/bit
Efs = 10*10^(-12);
Emp = 13*10^(-16);
EDA=5*10^(-9);


% Transmit Amplifier Types %
Eamp=100*10^(-12); % units in Joules/bit/m^2 (amount of energy spent by the amplifier to transmit the bits)


%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%

% Threshold Distance
do = sqrt(Efs/Emp);

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
    nodes(i).dist_from_origin=sqrt((nodes(i).x-sinkx)^2+(nodes(i).y-sinky)^2);
    nodes(i).dst_leader = 0;
    nodes(i).dst_leader = 0;
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


%G = graph();




operating_nodes = NUM_NODES;
dead_nodes=0;
round = 0;
while dead_nodes < NUM_NODES
    
    %M_S
    flag=0;
    viscircles([50,50],radius);
    if ((mod(round,4)==0))
        ms_Po.x = sinkx+radius;
        ms_Po.y = sinky;
    end
    if(mod(round,4)==1)
        ms_Po.x = sinkx;
        ms_Po.y = sinky+radius;
    end
    if(mod(round,4)==2)
        ms_Po.x = sinkx-radius;
        ms_Po.y = sinky;
    end
    if(mod(round,4)==3)
        ms_Po.x = sinkx;
        ms_Po.y = sinky-radius;
    end
    plot(ms_Po.x,ms_Po.y,'o','Linewidth',3);
    CLheads=0;
    for i = 1:NUM_NODES
        nodes(i).dst_ms = sqrt((nodes(i).x-ms_Po.x)^2+(nodes(i).y-ms_Po.y)^2);
    end
    %calculating weight
    for i = 1: k
        counter = 1;
        for j = 1 : NUM_NODES
            if(nodes(j).rleft == 0)
                nodes(j).weight = (nodes(j).battery)^2/nodes(j).dst_ms;
            
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
%         CH_s_dcir(i) = nodes(CH_s(i).id).dcir;   
        CLheads=CLheads+1;	% sum of cluster heads that have been elected 
        CL(CLheads).x=nodes(CH_s(i).id).x; % X-axis coordinates of elected cluster head
        CL(CLheads).y=nodes(CH_s(i).id).y; % Y-axis coordinates of elected cluster head
        CL(CLheads).id=CH_s(i).id; % Assigns the node ID of the newly elected cluster head to an array
        CL(CLheads).route = [];
        CL(CLheads).path = 0;
        CL(CLheads).dcir=sqrt((CL(CLheads).x-ms_Po.x)^2+(CL(CLheads).y-ms_Po.y)^2);
        DCIR(CLheads) = CL(CLheads).dcir;
        
        for j = 1:NUM_NODES
            if (nodes(j).cluster == i)
                nodes(j).dst_from_CH = sqrt((nodes(j).x-nodes(cluster(ID,i)).x)^2 + (nodes(j).y-nodes(cluster(ID,i)).y)^2);%Distance of node from its CH
            end
        end
    end
    % Fixing the size of "CL" array %
    CL=CL(1:CLheads);
    
    % Fixing the size of "DCIR" array
    DCIR = DCIR(1:CLheads);    
    [least_dcir,least_CH_id]=min(DCIR);
    
%     [leader_CH, ID_leader] = min(CH_s_dcir);
%     for i = 1:k
%         nodes(CH_s(i).id).dst_leader = sqrt((nodes(CH_s(i).id).x-nodes(CH_s(ID_leader).id).x)^2 + (nodes(CH_s(i).id).y-nodes(CH_s(ID_leader).id).y)^2 );
%     end
    


    
    

    %Graph
%     for i = 1:k
%         if(findnode(G,CH_s(i).id)==0 )
%             G=addnode(G,CH_s(i).id);
%         end
%         for j = i:k
%             if (i~=j)
%                 if((nodes(CH_s(j).id).dst_leader <= threshold)&&(findnode(G,CH_s(ID_leader).id)~=0)&&(findedge(G,CH_s(j).id,CH_s(ID_leader).id)==0))
%                     G = addedge(G, CH_s(ID_leader).id, CH_s(i).id, nodes(CH_s(j).id).dst_leader);
%                 elseif(findnode(G, CH_s(i).id)==0)
%                     G = addedge(G, CH_s(ID_leader).id, CH_s(i).id, nodes(CH_s(j).id).dst_from_CH);
%                 end
% 
%                 node_dist = ((nodes(i).x - nodes(j).x)^2 + (nodes(i).y - nodes(j).y)^2);
%                 if node_dist < threshold
%                     G = addedge(G, CH_s(i).id, CH_s(j).id, node_dist);
%                 end
%             end
%         end
%     end
    
%     for i = 1:k
%         nodes(CH_s(i).id).route = shortestpath(G, CH_s(i).id, CH_s(ID_leader).id);
%     end
    
    %%Energy Consumption during transmission and reception 
    for i = 1:NUM_NODES
        if (nodes(i).cond == 1)
            if (nodes(i).role == 0)
%                 nodes(i).battery = nodes(i).battery - (Eelec*packet_length + Eamp*packet_length*nodes(i).dst_from_CH^2);
              
                    
                if (nodes(i).dst_from_CH < do)
                    nodes(i).battery = nodes(i).battery - (packet_length*Eelec + Efs*packet_length*(nodes(i).dst_from_CH)^2);
                end
                if (nodes(i).dst_from_CH >= do)
                    nodes(i).battery = nodes(i).battery - (packet_length*Eelec + Emp*packet_length*(nodes(i).dst_from_CH)^4);
                end
            end
            if (nodes(i).role == 1)
%                nodes(i).battery = nodes(i).battery - (Eelec+EDA)*packet_length;
                nodes(i).battery = nodes(i).battery - packet_length*Eelec;
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
    
    % distance from CH to Leader CH
    for i = 1:CLheads
        CL(i).dist=sqrt((CL(least_CH_id).x-CL(i).x)^2 + (CL(least_CH_id).y-CL(i).y)^2);
    end
    
    % distance from CH to next nearest neighbour
    nearest_neighbour=zeros;
    for i=1:CLheads
        for j=i:CLheads
            if(i~=j)
                nearest_neighbour(i,j)=sqrt((CL(i).x-CL(j).x)^2 + (CL(i).y-CL(j).y)^2);
                nearest_neighbour(j,i)=sqrt((CL(i).x-CL(j).x)^2 + (CL(i).y-CL(j).y)^2); 
            end
            if(i==j)
                nearest_neighbour(i,j)=inf;
            end
        end
    end
    
    % Greedy Algorithm
    for i=1:CLheads
        if(i~=least_CH_id)
            [neigh_CHs_dis,neigh_CHs_id] = find_neigh_CHs(i,nearest_neighbour,CLheads);
            for j=1:CLheads
                if((CL(neigh_CHs_id(j)).dist<CL(i).dist)&&(CL(neigh_CHs_id(j)).path<2)&&(neigh_CHs_id(j)~=least_CH_id))
                            CL(i).path=CL(i).path+1;
                            CL(neigh_CHs_id(j)).path=CL(neigh_CHs_id(j)).path+1;
                            CL(i).route(length(CL(i).route)+1) = neigh_CHs_id(j); 
                            break;
                elseif((CL(neigh_CHs_id(j)).dist<CL(i).dist)&&(neigh_CHs_id(j)==least_CH_id))
                        CL(i).path=CL(i).path+1;
                        CL(neigh_CHs_id(j)).path=CL(neigh_CHs_id(j)).path+1;
                        CL(i).route(length(CL(i).route)+1) = neigh_CHs_id(j); 
                        break;
                end
            end
            else
                CL(i).route = i;
        end
    end 
    
    % Updating distance of CHs to closest CH
    for i=1:CLheads
        CL(i).dist=sqrt((CL(CL(i).route).x-CL(i).x)^2 + (CL(CL(i).route).y-CL(i).y)^2); 
    end
    % Energy reception of CHs due to Greedy Algorithm

    for i=1:CLheads
      if((nodes(CL(CL(i).route).id).battery>0)&&(nodes(CL(CL(i).route).id).cond==1)&&(nodes(CL(CL(i).route).id).role==1)&&(i~=least_CH_id))
          ERx=(Eelec)*packet_length;
%           energy=energy+ERx;
          nodes(CL(CL(i).route).id).battery=nodes(CL(CL(i).route).id).battery - ERx;
          if nodes(CL(CL(i).route).id).battery<=0  % if cluster heads energy depletes with reception
              nodes(CL(CL(i).route).id).cond=0;
              dead_nodes=dead_nodes +1;
              operating_nodes= operating_nodes - 1;
          end  
      end
    end
    
    % Energy Dissipation for cluster head nodes %
   
   for i=1:CLheads  
     if (nodes(CL(i).id).cond==1)  && (nodes(CL(i).id).role==1 ) && (i~=least_CH_id)
         if (nodes(CL(i).id).battery)>0
            if(CL(i).dist<do)
            ETx= (Eelec)*packet_length + Efs*packet_length*CL(i).dist^2;
            nodes(CL(i).id).battery=nodes(CL(i).id).battery- ETx;
%             energy=energy+ETx;
            else
            ETx= (Eelec)*packet_length + Emp*packet_length*CL(i).dist^4;
            nodes(CL(i).id).battery=nodes(CL(i).id).battery- ETx;
%             energy=energy+ETx;                
            end
         end
         if  nodes(CL(i).id).battery<=0     % if cluster heads energy depletes with transmission
         dead_nodes=dead_nodes +1;
         operating_nodes= operating_nodes - 1;
         nodes(CL(i).id).cond=0;
         end
     
     elseif(nodes(CL(i).id).cond==1)  && (nodes(CL(i).id).role==1 ) && (i==least_CH_id)
         if (nodes(CL(i).id).battery)>0
            if(CL(i).dist<do)
            ETx= (Eelec)*packet_length + Efs*packet_length*CL(i).dcir^2;
            nodes(CL(i).id).battery=nodes(CL(i).id).battery- ETx-packet_length*Eelec;
%             energy=energy+ETx;                
            else
            ETx= (Eelec)*packet_length + Emp*packet_length*CL(i).dcir^4;
            nodes(CL(i).id).battery=nodes(CL(i).id).battery- ETx-packet_length*Eelec;
%             energy=energy+ETx;                  
            end

         end
         if  nodes(CL(i).id).battery<=0     % if cluster heads energy depletes with transmission
         dead_nodes=dead_nodes +1;
         operating_nodes= operating_nodes - 1;
         nodes(CL(i).id).cond=0;
%          nodes(CL(i).id).rop=rnd;
         end
     end
   end
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
            %nodes(i).dst_leader = 0;
        end
    end
%     break; 
end

%Plotting Simulation Results "Operating Nodes per Round" %
    figure(2)
    plot(1:round,op(1:round),'-r','Linewidth',2);
    %axis([0  20000    0  NUM_NODES]);
    title ({'Updated k_means'; 'Operating Nodes per Round';})
    xlabel 'Rounds ';
    ylabel 'Operational Nodes ';
    hold on;
figure(3);    
% plot(G);
hold on;

% Function for finding nearest neighbour_CHs and Distance between them
function[neigh_CHs_dis,neigh_CHs_id] = find_neigh_CHs(i,nearest_neighbour,CLheads)
    neigh_CHs_id = [];
    neigh_CHs_dis = nearest_neighbour(i,:);
    neigh_CHs_dis = sort(neigh_CHs_dis);
    for j =1:CLheads
        for k=1:CLheads
            if(neigh_CHs_dis(j)==nearest_neighbour(i,k))
                neigh_CHs_id(length(neigh_CHs_id)+1) = k;
            end
        end
    end
end


