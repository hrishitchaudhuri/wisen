clc; clear all; close all;
NUM_NODES = 100;
no_of_clusters = 5;
angle_sector = 2*pi/no_of_clusters;
radius_field = 100;
x0 = 0;
y0 = 0;
packet_length = 500;
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
        %%Calculating weights and electing CH's
        for i = 1:no_of_clusters
            counter = 1;
            for j = 1:NUM_NODES
                if(nodes(j).cluster==i)
                    nodes(j).role=0;
                    nodes(j).weight = nodes(j).battery^2/nodes(j).dist_origin;
                    clusters(counter, i)=nodes(j).id;
                    weights(counter, i)=nodes(j).weight;
                    counter = counter+1;
                end
                CH_s(i).no_of_nodes = counter;
            end
            [weight, ID] = max(weights(:,i));
            CH_s(i).id = clusters(ID,i) ;
            nodes(CH_s(i).id).role = 1;
        end

        %%Calculating distance of nodes from CH's
        for i = 1:no_of_clusters
            for j = 1:NUM_NODES
                if(nodes(j).cluster == i)
                    nodes(j).dist_CH = sqrt((nodes(j).x-nodes(CH_s(i).id).x)^2+(nodes(j).y-nodes(CH_s(i).id).y)^2);
                end
            end
        end
        %Energy dissipation
        
        for i = 1:NUM_NODES
            if(nodes(i).cond == 1)
                if (nodes(i).role == 0)
                    if (nodes(i).dist_CH < do)
                        nodes(i).battery = nodes(i).battery - (packet_length*Eelec + Efs*packet_length*(nodes(i).dist_CH)^2);
                    end
                    if (nodes(i).dist_CH >= do)
                        nodes(i).battery = nodes(i).battery - (packet_length*Eelec + Emp*packet_length*(nodes(i).dist_CH)^4);
                    end
                end
                if(nodes(i).role == 1)
                    %For reception
                    nodes(i).battery = nodes(i).battery - packet_length*(CH_s(nodes(i).cluster).no_of_nodes-1)*Eelec;
                    %For transmission to static sink
                    if (nodes(i).battery > 0)
                        if (nodes(i).dist_origin < do)
                            nodes(i).battery = nodes(i).battery - (CH_s(nodes(i).cluster).no_of_nodes*packet_length*Eelec + Efs*CH_s(nodes(i).cluster).no_of_nodes*packet_length*(nodes(i).dist_origin)^2);
                        end
                        if (nodes(i).dist_origin >= do)
                            nodes(i).battery = nodes(i).battery - (CH_s(nodes(i).cluster).no_of_nodes*packet_length*Eelec + Emp*CH_s(nodes(i).cluster).no_of_nodes*packet_length*(nodes(i).dist_origin)^4);
                        end
                    end
                end
            if(nodes(i).battery<=0)
                nodes(i).cond = 0;
%                 nodes(i).cluster = 0;
                nodes(i).rip = rounds;
                dead_nodes = dead_nodes + 1;
                operating_nodes = operating_nodes - 1;
                disp(nodes(i).id);
            end
            end
         end
        rounds = rounds+1;
        disp(rounds);
        op(rounds)=operating_nodes;
        
%         Resetting variables for next round 
%         clusters = zeros;
%         weights = zeros;
    end

    
figure(2)
plot(1:rounds,op(1:rounds),'-g','Linewidth',2);
% axis([0  6000    0  NUM_NODES]);

hold on;
    