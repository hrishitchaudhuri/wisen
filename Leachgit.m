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
%Optimal Election Probability of a node to become cluster head/ 
p=0.1;
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
%Values for Hetereogeneity 
%Percentage of nodes than are advanced 
% m=0.1;
%\alpha
a=1;
%%%sensor data%%%
%temperature measurement
tempi=50; %initial temperature
tempf=200; %final temperature
sv=0; %previously sensed value
cv=0; %current sensing value

%Threshold for transmiting data to the cluster head
h=100;                                 %%%%%%Hard Threshold H(t)
s=2;                                   %%%%%%Soft threshold  S(t)
%maximum number of rounds
rmax=15;
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
%Computation of do
do=sqrt(Efs/Emp);
%%%%%% Initialising the arrays %%%%%%%
PACKETS_TO_CH=zeros;
PACKETS_TO_BS=zeros;
% DEAD_N=zeros;
X=zeros;
Y=zeros;
sdata=zeros;
% C=zeros;
% STATISTICS=zeros;
% CLUSTERHS=zeros;
% CH_dis=zeros;
%%%%%%%%%%%%%%%%%
%Creation of the random Sensor Network
figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';
   
%     temp_rnd0=i;
%     Replaced i instead of temp_rand the if condition. 
    %Random Election of Normal Nodes 
%     if (i>=m*n+1) 
        S(i).E=Eo;
%         S(i).ENERGY=0;
%         plot(S(i).xd,S(i).yd,'red o');
%         hold on;
%     end
    %Random Election of Advanced Nodes
%    if (temp_rnd0<m*n+1)  
%         S(i).E=Eo*(1+a)
%         S(i).ENERGY=1;
%         plot(S(i).xd,S(i).yd,'+');
%         hold on;
%     end
end
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
% plot(S(n+1).xd,S(n+1).yd,'green x');
% hold on;
  
      
%First Iteration
figure(1);
%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;
countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;
  for r=1:1:rmax
%     r=0
    disp(r);
%     hold on;
    plot(S(n+1).xd,S(n+1).yd,'green x');
    hold on;
    
  %Operation for epoch
   if(mod(r,round(1/p))==0)
     for i=1:1:n
         S(i).G=0;
         S(i).cl=0; 
     end
   end
% hold off;
%Number of dead nodes
% dead=0;
%Number of dead Advanced Nodes
% dead_a=0;
%Number of dead Normal Nodes
dead_n=0;
%counter for bit transmitted to Bases Station and to Cluster Heads 
packets_TO_BS=0;
packets_TO_CH=0;
%counter for bit transmitted to Bases Station and to Cluster Heads 
%per round
PACKETS_TO_CH(r)=0;
PACKETS_TO_BS(r)=0;
figure(1);
for i=1:1:n
    %checking if there is a dead node 
%     if (S(i).E<=0)
%         plot(S(i).xd,S(i).yd,'red +');
%         dead=dead+1;
%         if(S(i).ENERGY==1)
%             dead_a=dead_a+1;
%         end
        %changed Energy to E
        if(S(i).E==0)
            dead_n=dead_n+1;
            plot(S(i).xd,S(i).yd,'red +');
        end
%         hold on;    
%     end
     if (S(i).E>0)
         S(i).type='N';
%         if (S(i).ENERGY==0)  
         plot(S(i).xd,S(i).yd,'red o');
          hold on;
     end
%         if (S(i).ENERGY==1)  
%         plot(S(i).xd,S(i).yd,'+');
%         end
%         hold on;
%     end
end
% plot(S(n+1).xd,S(n+1).yd,'x');
% STATISTICS(r+1).DEAD=dead;
% DEAD(r+1)=dead;
DEAD_N(r)=dead_n;
% DEAD_A(r+1)=dead_a;
%When the first node dies
if (dead_n==1)
    if(flag_first_dead==0)
        first_dead=r;
        flag_first_dead=1;
    end
end
countCHs=0;
cluster=1;
for i=1:1:n
   if(S(i).E>0)
   temp_rand=rand;     
   if ((S(i).G)<=0)
 %Selection of Cluster Heads/ 
 if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
            countCHs=countCHs+1;
%figure the logic for packets to BS
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r)=packets_TO_BS;
            
            S(i).type='C';
            S(i).G=round(1/p); 
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            plot(S(i).xd,S(i).yd,'k*');
            
            distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
            C(cluster).distance=distance;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            cluster=cluster+1;
            
            %Calculation of Energy dissipated for CH 
            if (distance>do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance^4)); 
            end
            if (distance<=do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*(distance^2)); 
            end
 end     
    
    end
  end 
end
%Delete this
cluster=cluster-1;
STATISTICS(r).CLUSTERHEADS=cluster;
CLUSTERHS(r)=cluster;
%Election of Associated Cluster Head for Normal Nodes
% min_dis ---> base_dis
% temp--->CH_dis

[base_dis,min_dis_cluster,energy_nodes,packets_TO_CH,S]=evaluate_CHs(S,C,cluster,packets_TO_CH);


% for i=1:1:n
%    if ( S(i).type=='N' && S(i).E>0 )
%      if(cluster-1>=1)
%        base_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
%        min_dis_cluster=1;
%         for CHcounter=1:1:cluster
%             CH_dis(CHcounter)=sqrt( (S(i).xd-C(CHcounter).xd)^2 + (S(i).yd-C(CHcounter).yd)^2);
%         end
%         [CH_dis,min_dis_cluster]=min(CH_dis);
%         S(i).cluster = min_dis_cluster;
%         if (CH_dis<base_dis)
%             base_dis=CH_dis;
%         end
        
     
%        for c=1:1:cluster-1
%            CH_dis=min(base_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
%            if ( CH_dis<base_dis )
%                base_dis=CH_dis;
%                min_dis_cluster=c;
%            end
%        end
       
       %Energy dissipated by sending packets to Cluster Head by normal
       %nodes 
            %if (cv >= h)
            %test = cv-sv;
            %if (test >= s)
       
%             base_dis;
%             if (base_dis>do)
%                 S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( base_dis^4)); 
%             end
%             if (CH_dis<=do)
%                 S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( base_dis^2)); 
%             end
        %Energy dissipated by cluster heads
%         if(base_dis>0)
%           S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
%          packets_TO_CH = (packets_TO_CH)+1; 
%         end
%        S(i).base_dis=base_dis;
%        S(i).cluster=min_dis_cluster;
           
%      end
            %end
            %end
            %    end
%    energy_nodes(i)= S(i).E;
% end
 
% countCHs;
rcountCHs=rcountCHs+countCHs;
%     for i=1:1:n
%         if(S(i).type=='C')
%             S(i).G = S(i).G-1;
    %    S(i).cl=0; 
%         end
%     end
PACKETS_TO_CH(r)=packets_TO_CH;
  for i=1:1:n
      if(S(i).type=='N')
      if(S(i).cluster==1)
          plot(S(i).xd,S(i).yd,'red o');
      end
      if(S(i).cluster==2)
          plot(S(i).xd,S(i).yd,'green o');
      end
      if(S(i).cluster==3)
          plot(S(i).xd,S(i).yd,'blue o');
      end
      if(S(i).cluster==4)
          plot(S(i).xd,S(i).yd,'yellow o');
      end
       if(S(i).cluster==5)
           plot(S(i).xd,S(i).yd,'red s');
       end
      if(S(i).cluster==5)
          plot(S(i).xd,S(i).yd,'magenta o');
      end
      if(S(i).cluster==6)
          plot(S(i).xd,S(i).yd,'black o');
      end
      if(S(i).cluster==7)
          plot(S(i).xd,S(i).yd,'cyan o');
      end
       if(S(i).cluster==8)
          plot(S(i).xd,S(i).yd,'blue d');
       end
       if(S(i).cluster==9)
          plot(S(i).xd,S(i).yd,'magenta ^');
       end
       if(S(i).cluster==10)
          plot(S(i).xd,S(i).yd,'black <');
       end
       if(S(i).cluster==11)
          plot(S(i).xd,S(i).yd,'cyan v');
       end
       if(S(i).cluster==12)
          plot(S(i).xd,S(i).yd,'red ^');
      end
  end
  end
  hold off;
  end
%%%%%%%%%%%
%%% Functions %%%%%

    






%%%%%%%%%%%%
% Path for mobile sink
hold on
%set axis limits ahead of time
axis([0 100 0 100]); 
%sorting
for i=1:countCHs
    for j=i+1:countCHs
        if (X(i)>X(j))
            a =  X(i);
            b =  Y(i);
            X(i) = X(j);
            Y(i) = Y(j);
            X(j) = a;
            Y(j) = b;
        end
    end
end
% plot empty line objects (with NaN values)
line1 = plot(X, nan(size(X)),'-','Color','r');
marker1 = plot(nan,nan,'*','Color','r');
% Path for joining CHs
for k=1:countCHs
    %marker plots
    marker1.XData = X(k);
    marker1.YData = Y(k);
    
    %Line plots
    line1.YData(k) = Y(k);
    
    drawnow(); 
    pause(0.5);
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   STATISTICS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     %
%  DEAD  : a rmax x 1 array of number of dead nodes/round                             %
%  DEAD_A : a rmax x 1 array of number of dead Advanced nodes/round                   %
%  DEAD_N : a rmax x 1 array of number of dead Normal nodes/round                     %
%  CLUSTERHS : a rmax x 1 array of number of Cluster Heads/round                      %
%  PACKETS_TO_BS : a rmax x 1 array of number packets send to Base Station/round      %
%  PACKETS_TO_CH : a rmax x 1 array of number of packets send to ClusterHeads/round   %
%  first_dead: the round where the first node died                                    %
%                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=1:rmax;
i=1:n;
figure (2);
plot(r,PACKETS_TO_CH(r));
xlabel('Rounds');
ylabel('Number of packets sent to the CHs');
title('LEACH');
figure(3);
plot(r,PACKETS_TO_BS(r));
xlabel('Rounds');
ylabel('Number of packets sent by the CHs to BS');
title('LEACH');
figure(4);
plot(r,CLUSTERHS(r));
xlabel('Rounds');
ylabel('Number of cluster heads');
title('LEACH');
figure(5);
plot(i,energy_nodes);
xlim([1 50]);
ylim([0.2 0.9]);
xlabel('Nodes');
ylabel('Energy');
title('LEACH');


%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%

% Evaluates CHs and calculates energy to the CH
function [base_dis,min_dis_cluster,energy_nodes,packets_TO_CH,S] = evaluate_CHs(S,C,cluster,packets_TO_CH)
%     Transmit Amplifier types 
    Efs=10*0.000000000001;
    Emp=0.0013*0.000000000001;
%     Computation of do
    do=sqrt(Efs/Emp);
%     %Eelec=Etx=Erx 
    ETX=50*0.000000001;
    ERX=50*0.000000001;
%     Data Aggregation Energy/ 
    EDA=5*0.000000001;
    n=50;
    for i=1:1:n
    if ( S(i).type=='N' && S(i).E>0 )
     if(cluster>=1)
           base_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
           min_dis_cluster=1;
           for CHcounter=1:1:cluster
            CH_dis(CHcounter)=sqrt( (S(i).xd-C(CHcounter).xd)^2 + (S(i).yd-C(CHcounter).yd)^2);
           end
           [CH_dis,min_dis_cluster]=min(CH_dis);
%         S(i).cluster = min_dis_cluster;
           if (CH_dis<base_dis)
            base_dis=CH_dis;
            if (base_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( base_dis^4)); 
            end
            if (CH_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( base_dis^2)); 
            end
        %Energy dissipated by cluster heads
            if(base_dis>0)
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
                packets_TO_CH = (packets_TO_CH)+1; 
            end
           end
     S(i).base_dis=base_dis;
     S(i).cluster=min_dis_cluster;
     end
   energy_nodes(i)= S(i).E;
    end
     end
end

%%%%Function to get sensor data (temperature) for each node during a
%%%%particular round%%%%%%%
function [sdata,cv,sv] = snodevalues(tempi,tempf,h,s)
for i=1:1:n
        cv = tempi + (tempf-tempi).*rand(1,1); 
        sdata(i) = cv;
        disp(sdata(i));
    end
end