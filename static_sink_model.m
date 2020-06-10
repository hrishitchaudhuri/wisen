clc;clear;close;
%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;
%x and y Coordinates of the Sink  
sink.x=0.5*xm;
sink.y=0.5*ym;
%Number of Nodes in the field  
n=50;
%Energy Model (all values in Joules)
%Initial Energy 
Eo=0.5;
%Eelec=Etx=Erx 
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types 
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
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
X=zeros;
Y=zeros;
figure(1);
hold on;
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    plot(S(i).xd, S(i).yd, "green o");
    S(i).E=Eo;
end
plot(sink.x, sink.y, "red ^");
hold off;
%Energy dissipated by the nodes
for r=1:25
    for i=1:n
        distance=sqrt((S(i).xd-sink.x)^2 + (S(i).yd-sink.y)^2 );
        if (distance>do)
                    S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*(distance^4)); 
        end
        %% 
        % DESCRIPTIVE TEXT
        if (distance<=do)
                    S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*(distance^2)); 
        end
        energy_nodes(i)= S(i).E;
        disp(S(i).E);
    end
end
figure (2);
i=1:n;
plot(i, energy_nodes);
xlim([1 50]);
ylim([0.45 0.55]);
xlabel('Nodes');
ylabel('Energy');
title('Static sink');

    