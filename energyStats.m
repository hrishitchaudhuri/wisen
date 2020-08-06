clc; clear all; close all;
NUM_NODES = 100;
no_of_clusters = 5;
angle_sector = 2*pi/no_of_clusters;
radius_field = 100;
x0 = 0;
y0 = 0;
packet_length = 500;
ad_length = 10;

radius_ms = 25;
%%%%Energy parameters
Eo = 0.5;
Eelec=50*10^(-9); % units in Joules/bit
Efs = 10*10^(-12);
Emp = 13*10^(-16);
EDA=5*10^(-9);
do = sqrt(Efs/Emp);


figure(1);
viscircles([x0,y0],radius_field);
hold on
for i = 1:NUM_NODES
    t = 2*pi*rand(1,1);
    r = radius_field*sqrt(rand(1,1));
    S(i, 1) = x0 + r*cos(t);
    S(i, 2) = y0 + r*sin(t);
    plot(S(i, 1), S(i, 2), 'red .');
    title 'Wireless Sensor Network';
    xlabel '(m)';
    ylabel '(m)';
    hold on;
end

a2 = 0;
for i = 1:no_of_clusters
    a1 = a2;  % A random direction
    a2 = a1 + angle_sector;
    t = linspace(a1,a2);
    x = x0 + radius_field*cos(t);
    y = y0 + radius_field*sin(t);
    plot([x0,x,x0],[y0,y,y0],'k -')
    axis equal
    hold on;
end

for i = 1:NUM_NODES
    nodes(i).id=i;	% sensor's ID number
    nodes(i).x=S(i, 1);	% X-axis coordinates of sensor node
    nodes(i).y=S(i, 2);	% Y-axis coordinates of sensor node
    nodes(i).battery=Eo;     % nodes energy levels (initially set to be equal to "Eo"
    nodes(i).role=0;   % node acts as normal if the value is '0', if elected as a cluster head it  gets the value '1' (initially all nodes are normal)
    nodes(i).cluster=0;	% the cluster which a node belongs to
    nodes(i).cond=1;
    nodes(i).dist_origin = sqrt((nodes(i).x-x0)^2+(nodes(i).y-y0)^2);    
    nodes(i).role=0;
    nodes(i).cond=1;
    nodes(i).dist_CH=0;
end

for i = 1:NUM_NODES
    if (nodes(i).x>=0&&nodes(i).y>=0)
        nodes(i).cluster = ceil((atan(nodes(i).y/nodes(i).x))/angle_sector);
    elseif(nodes(i).x<0)
        nodes(i).cluster = ceil((atan(nodes(i).y/nodes(i).x)+pi)/angle_sector);
    elseif(nodes(i).x>=0&&nodes(i).y<0)
        nodes(i).cluster = ceil((atan(nodes(i).y/nodes(i).x)+2*pi)/angle_sector);
    end
end






    dead_nodes = 0;
    rounds = 0;
    operating_nodes = NUM_NODES;
    while dead_nodes<NUM_NODES
        flag=0;
        viscircles([x0,y0],radius_ms);
        if ((mod(rounds,4)==0))
            ms_Po.x = x0+radius_ms;
            ms_Po.y = y0;
        end
        if(mod(rounds,4)==1)
            ms_Po.x = x0;
            ms_Po.y = y0+radius_ms;
        end
        if(mod(rounds,4)==2)
            ms_Po.x = x0-radius_ms;
            ms_Po.y = y0;
        end
        if(mod(rounds,4)==3)
            ms_Po.x = x0;
            ms_Po.y = x0-radius_ms;
        end
        plot(ms_Po.x,ms_Po.y,'o','Linewidth',3);
        for i = 1:NUM_NODES
            nodes(i).dist_ms = sqrt((nodes(i).x-ms_Po.x)^2+(nodes(i).y-ms_Po.y)^2);
        end
        %%Calculating weights and electing CH's
        for i = 1:no_of_clusters
            counter = 1;
            CH_s(i).path = 0;
            CH_s(i).route = [];
            for j = 1:NUM_NODES
                if(nodes(j).cluster==i && nodes(j).cond == 1)
                    nodes(j).role=0;
                    nodes(j).weight = nodes(j).battery^2/nodes(j).dist_ms;
                    clusters(counter, i)=nodes(j).id;
                    weights(counter, i)=nodes(j).weight;
                    CH_s(i).no_of_nodes = counter;
                    counter = counter+1;
                end
                
            end
            if(CH_s(i).no_of_nodes ~= 0 )
                [weight, ID] = max(weights(:,i));
                CH_s(i).id = clusters(ID,i);
                if (nodes(clusters(ID,i)).cond ==1)
                   nodes(CH_s(i).id).role = 1; 
                   CH_each_round(i, rounds+1) = CH_s(i).id;
                   CH_s(i).x = nodes(CH_s(i).id).x;
                   CH_s(i).y = nodes(CH_s(i).id).y;
                   CH_s(i).dist_ms = nodes(CH_s(i).id).dist_ms;
                   dist_MS(i) = CH_s(i).dist_ms;
                end
            else
                CH_s(i).id = 0;
                CH_s(i).dist_ms = inf;
                CH_s(i).dist = inf;
                dist_MS(i) = CH_s(i).dist_ms;
                CH_s(i).path = inf;
                CH_s(i).route = inf;                
                CH_each_round(i, rounds+1) = CH_s(i).id;
            end
            
        end
        %%Leader CH election
        [leader_dist_MS, leader_CH_ID] = min(dist_MS);
        
        
        %%Calculating distance of nodes from CH's
        for i = 1:no_of_clusters
            for j = 1:NUM_NODES 
                if(nodes(j).cluster == i && nodes(j).cond == 1 )
                    nodes(j).dist_CH = sqrt((nodes(j).x-nodes(CH_s(i).id).x)^2+(nodes(j).y-nodes(CH_s(i).id).y)^2);
                    if (nodes(i).cond == 1 && nodes(i).role == 1)
                       if (nodes(j).dist_CH < do)
                            nodes(i).battery = nodes(i).battery - (ad_length*Eelec + Efs*ad_length*(nodes(j).dist_CH)^2);
                        end
                        if (nodes(i).dist_CH >= do)
                            nodes(i).battery = nodes(i).battery - (ad_length*Eelec + Emp*ad_length*(nodes(j).dist_CH)^4);
                        end 
                         if(nodes(i).battery<=10^-3)
                            nodes(i).cond = 0;
                            nodes(i).role = 0;
                            
            %                 nodes(i).cluster = 0;
                            nodes(i).rip = rounds;
                            dead_nodes = dead_nodes + 1;
                            operating_nodes = operating_nodes - 1;
                            disp(nodes(i).id);
                            CH_s(nodes(i).cluster).no_of_nodes = CH_s(nodes(i).cluster).no_of_nodes -1 ;
    
                        end
                    end 
                end 
            end  
        end
        %Energy dissipation
        
        for i = 1:NUM_NODES
            if(nodes(i).cond == 1)
                if (nodes(i).role == 0 )
                    if (nodes(i).dist_CH < do)
                        nodes(i).battery = nodes(i).battery - (packet_length*Eelec + Efs*packet_length*(nodes(i).dist_CH)^2);
                    end
                    if (nodes(i).dist_CH >= do)
                        nodes(i).battery = nodes(i).battery - (packet_length*Eelec + Emp*packet_length*(nodes(i).dist_CH)^4);
                    end
                    nodes(i).battery = nodes(i).battery - ad_length*Eelec;
                end
                if(nodes(i).role == 1)
                    %For reception
                    if(CH_s(nodes(i).cluster).no_of_nodes>1)
                        nodes(i).battery = nodes(i).battery - packet_length*(CH_s(nodes(i).cluster).no_of_nodes-1)*Eelec;
                    else
                        nodes(i).battery = nodes(i).battery - packet_length*(CH_s(nodes(i).cluster).no_of_nodes)*Eelec;
                    end
                    %For transmission to static sink
%                     if (nodes(i).battery > 0)
%                         if (nodes(i).dist_origin < do)
%                             nodes(i).battery = nodes(i).battery - (CH_s(nodes(i).cluster).no_of_nodes*packet_length*Eelec + Efs*CH_s(nodes(i).cluster).no_of_nodes*packet_length*(nodes(i).dist_origin)^2);
%                         end
%                         if (nodes(i).dist_origin >= do)
%                             nodes(i).battery = nodes(i).battery - (CH_s(nodes(i).cluster).no_of_nodes*packet_length*Eelec + Emp*CH_s(nodes(i).cluster).no_of_nodes*packet_length*(nodes(i).dist_origin)^4);
%                         end
%                     end
                end
            if(nodes(i).battery<=10^-3)
                nodes(i).cond = 0;
                nodes(i).role = 0;
                
%                 nodes(i).cluster = 0;
                nodes(i).rip = rounds;
                dead_nodes = dead_nodes + 1;
                operating_nodes = operating_nodes - 1;
%                 disp(nodes(i).id);
                CH_s(nodes(i).cluster).no_of_nodes = CH_s(nodes(i).cluster).no_of_nodes -1 ;

            end
            end
        end
        
     % distance from CH to Leader CH
    for i = 1:no_of_clusters
        if(CH_s(i).id>0)
            CH_s(i).dist=sqrt((CH_s(leader_CH_ID).x-CH_s(i).x)^2 + (CH_s(leader_CH_ID).y-CH_s(i).y)^2);
        end
    end
    
    % distance from CH to next nearest neighbour
    nearest_neighbour=zeros;
    for i=1:no_of_clusters
        for j=i:no_of_clusters
            if(i~=j&&CH_s(i).id>0&&CH_s(j).id>0)
                nearest_neighbour(i,j)=sqrt((CH_s(i).x-CH_s(j).x)^2 + (CH_s(i).y-CH_s(j).y)^2);
                nearest_neighbour(j,i)=sqrt((CH_s(i).x-CH_s(j).x)^2 + (CH_s(i).y-CH_s(j).y)^2); 
            end
            if(i==j||CH_s(i).id==0||CH_s(j).id==0)
                nearest_neighbour(i,j)=inf;
                nearest_neighbour(j,i)=inf;
            end
        end
    end
    
    % Greedy Algorithm
    for i=1:no_of_clusters
        if(i~=leader_CH_ID&&CH_s(i).id>0)
            [neigh_CHs_dis,neigh_CHs_id] = find_neigh_CHs(i,nearest_neighbour,no_of_clusters,CH_s);
            for j=1:no_of_clusters
                if(CH_s(neigh_CHs_id(j)).id>0)
                if((CH_s(neigh_CHs_id(j)).dist< CH_s(i).dist)&&(CH_s(neigh_CHs_id(j)).path<2)&&(neigh_CHs_id(j)~=leader_CH_ID))
                        CH_s(i).path=CH_s(i).path+1;
                            CH_s(neigh_CHs_id(j)).path=CH_s(neigh_CHs_id(j)).path+1;
                        CH_s(i).route(length(CH_s(i).route)+1) = neigh_CHs_id(j); 
                        break;
                elseif((CH_s(neigh_CHs_id(j)).dist<CH_s(i).dist)&&(neigh_CHs_id(j)==leader_CH_ID))
                        CH_s(i).path=CH_s(i).path+1;
                        CH_s(neigh_CHs_id(j)).path=CH_s(neigh_CHs_id(j)).path+1;
                        CH_s(i).route(length(CH_s(i).route)+1) = neigh_CHs_id(j); 
                        break;
                end
                end
            end
            
        elseif(i==leader_CH_ID)
                CH_s(i).route = i;
        end
     end
    
    % Updating distance of CHs to closest CH
    for i=1:no_of_clusters
        if(CH_s(i).id>0)
            CH_s(i).dist=sqrt((CH_s(CH_s(i).route).x-CH_s(i).x)^2 + (CH_s(CH_s(i).route).y-CH_s(i).y)^2); 
        end
    end
        
    for i=1:1:no_of_clusters
        j = i;
        transmit_nodes = CH_s(i).no_of_nodes;
        
        if j == leader_CH_ID && CH_s(i).id > 0
            nodes(CH_s(j).id).battery = nodes(CH_s(j).id).battery - packet_length*transmit_nodes*Eelec;
        end
                
        while j ~= leader_CH_ID && CH_s(j).id > 0 && j ~= Inf
            if CH_s(j).dist < do
                % fill stats 1
                nodes(CH_s(j).id).battery = nodes(CH_s(j).id).battery - (packet_length*transmit_nodes*Eelec + Efs*packet_length*transmit_nodes*(CH_s(j).dist)^2);
            else
                % fill stats 2
                nodes(CH_s(j).id).battery = nodes(CH_s(j).id).battery - (packet_length*transmit_nodes*Eelec + Emp*packet_length*transmit_nodes*(CH_s(j).dist)^4);
            end

            nodes(CH_s(CH_s(j).route).id).battery = nodes(CH_s(CH_s(j).route).id).battery - packet_length*transmit_nodes*Eelec;
            j = CH_s(j).route;
        end
    end
        
        
        
        
        
        rounds = rounds+1;
        disp(rounds);
        op(rounds)=operating_nodes;
        
        

%         Resetting variables for next round 
        clusters = zeros;
        weights = zeros;
        dist_MS = zeros;
       
    end

    
figure(2)
plot(1:rounds,op(1:rounds),'-g','Linewidth',2);
% axis([0  6000    0  NUM_NODES]);

hold on;


function[neigh_CHs_dis,neigh_CHs_id] = find_neigh_CHs(i,nearest_neighbour,no_of_clusters,CH_s)
    neigh_CHs_id = [];
    neigh_CHs_dis = nearest_neighbour(i,:);
    neigh_CHs_dis = sort(neigh_CHs_dis);
    for j =1:no_of_clusters
        for k=1:no_of_clusters
            if(neigh_CHs_dis(j)==nearest_neighbour(i,k))
                neigh_CHs_id(length(neigh_CHs_id)+1) = k;
            end
        end
    end
end