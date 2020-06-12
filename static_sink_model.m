clc;clear;close;
%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Field Dimensions - x and y maximum (in meters)
xm=100;
% %% 
%% 
ym=100;
%x and y Coordinates of the Sink  
sink.x=0.5*xm;
sink.y=0.5*ym;
%Number of Nodes in the field  
NUM_NODES=50;
NUM_ROUNDS = 25
%Energy Model (all values in Joules)
%Initial Energy 
Eo=0.5;
%Eelec=Etx=Erx 
ETX=50*0.0000001;
ERX=50*0.0000001;
%Transmit Amplifier types 
Efs=10*0.0000000001;
Emp=0.0013*0.0000000001;
%Data Aggregation Energy/ 
EDA=5*0.000000001;
%\alpha
%% 
a=1;
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
%Computation of do
do=sqrt(Efs/Emp);
%%
% 
%   for x = 1:10
%       disp(x)
%   end
% 
%% 
%% 
X=zeros;
Y=zeros;

[nodes,pos_nodes] = initialise(NUM_NODES);
plot(sink.x, sink.y, "red ^");
hold off;
%Energy dissipated by the nodes
[energy_nodes, energy_nodes_rounds] = energy_parameters(NUM_NODES, nodes, sink, NUM_ROUNDS);

node_num = input("Enter node number for viewing its energy statics");
figure (2);
i=1:NUM_NODES;
plot(i, energy_nodes_rounds(node_num));
%xlim([1 50]);
%ylim([0.45 0.55]);
xlabel('Rounds');
ylabel('Energy');
title('Static sink');



%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%
function[nodes,pos_nodes] = initialise(no_of_nodes)
    figure(1);
    hold on;
     for i=1:1:no_of_nodes
        nodes(i) = Sensor_Node;
        nodes(i)= Add_Positions(nodes(i));
        pos_nodes(i,1)=nodes(i).xd;
        pos_nodes(i,2)=nodes(i).yd;
        plot(nodes(i).xd,nodes(i).yd,'green o');
     end
end

function [energy_nodes, energy_nodes_rounds] = energy_parameters(no_of_nodes, NODES, sink, NUM_ROUNDS)
    Eo=0.5;
    %Eelec=Etx=Erx 
    ETX=50*0.000000001;
    ERX=50*0.00000001;
    %Transmit Amplifier types 
    Efs=10*0.0000000001;
    Emp=0.0013*0.0000000001;
    %Data Aggregation Energy/ 
    EDA=5*0.000000001;
    do=sqrt(Efs/Emp);
    for r=1:NUM_ROUNDS
        for i=1:no_of_nodes
            distance=sqrt((NODES(i).xd-sink.x)^2 + (NODES(i).yd-sink.y)^2 );
            if (distance>do)
                        NODES(i).E=NODES(i).E- ( ETX*(4096) + Emp*4096*(distance^4)); 
            end
            %% 
            % DESCRIPTIVE TEXT
            if (distance<=do)
                        NODES(i).E=NODES(i).E- ( ETX*(4096) + Efs*4096*(distance^2)); 
            end
            energy_nodes(i)= NODES(i).E;
            energy_nodes_rounds(i)= i, NODES(i).E;
            disp(NODES(i).E);
        end
    end
end


    
