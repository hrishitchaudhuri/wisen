clc;clear;close;

%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%

% Field Dimensions
x_max = 100;
y_max = 100;

% Sink Coordinates
sink.x = 0.5 * x_max;
sink.y = 0.5 * y_max;

% Number of Nodes
NUM_NODES = 5;

% Number of Clusters
global k;
k = 1;

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
rounds = 10;

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


[nodes] = assign_CHs(nodes,NUM_NODES,C,idx);

wake_clust = 1;
dead_nodes = 0;
plot_idx = 1;
round_stat = zeros;
dead_stats = zeros;

[nodes,dead_nodes,round_stat,dead_stats,wake_clust,distance] = energy_model(nodes,NUM_NODES,wake_clust,plot_idx,dead_nodes);

% [matr ix] = distance_between_nodes(nodes,NUM_NODES);
% disp(matrix);
% DISTANCE BASED ROUTING
% [go,trace] = distance_based_routing(nodes,matrix,NUM_NODES);

[nodes,source_nodes,G] = route_taken(nodes,NUM_NODES,distance);

figure(2);
plot(round_stat, dead_stats);
hold on;

figure(3);
plot(G);
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
    global k;
    for i = 1:1:NUM_NODES
    nodes(i).cluster = (idx(i)-1);
    nodes(i).sink_x = C(idx(i), 1);
    nodes(i).sink_y = C(idx(i), 2);
    end
    
    for i=NUM_NODES+1:1:NUM_NODES+k
        nodes(i).x = C(i-NUM_NODES,1);
        nodes(i).y = C(i-NUM_NODES,2);
        nodes(i).distance=0;
    end
    
end


function[nodes,dead_nodes,round_stat,dead_stats,wake_clust,distance] = energy_model(nodes,NUM_NODES,wake_clust,plot_idx,dead_nodes)
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
            
            distance(j)=sqrt((nodes(j).x - nodes(j).sink_x)^2 + (nodes(j).y - nodes(j).sink_y)^2);
            
            nodes(j).distance=distance(j);
            if (distance(j)>do)
                nodes(j).battery = nodes(j).battery - ((ETX)*(4096) + Emp*4096*(distance(j)^4)); 
            end
        
            if (distance(j)<=do)
               nodes(j).battery = nodes(j).battery - ((ETX)*(4096) + Efs*4096*(distance(j)^2)); 
            end
        
            if (nodes(j).battery <= 0)
                dead_nodes = dead_nodes + 1;
                nodes(j).state = 0;
                plot(nodes(j).x, nodes(j).y, 'black x');
            end
            
        end
    end
    
%      pause(0.5);
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

function[distance_nodes] = distance_between_nodes(nodes,i,j)
    distance_nodes = sqrt((nodes(i).x-nodes(j).x)^2+(nodes(i).y-nodes(j).y)^2);
end


% function[matrix] = distance_between_nodes(nodes,NUM_NODES)
%     global k;
%     for i=1:1:NUM_NODES+k
%         for j=1:1:NUM_NODES+k
%             matrix(i,j) = (nodes(i).x-nodes(j).x)^2+(nodes(i).y-nodes(j).y)^2;
%         end
%     end
%     
% end

%%%%%%%%%%%%%%%%%%%%%%%

function[nodes,source_nodes,G] = route_taken(nodes,NUM_NODES,distance)
     source_nodes = zeros;
     G=digraph();
     threshold = 12; 
     for i=1:1:NUM_NODES
        for j=1:1:NUM_NODES
            if((nodes(i).cluster==nodes(j).cluster)&&(i~=j))
                source_nodes(i,j) = j;
                [distance_nodes] = distance_between_nodes(nodes,i,j);
                if((distance_nodes<=threshold)&&(nodes(j).state==1))
                    G=addedge(G,i,j,distance_nodes);
                    if(nodes(j).distance<=threshold)
                        G=addedge(G,j,(nodes(j).cluster+NUM_NODES+1),nodes(j).distance);
                    else
                        G=addedge(G,j,(nodes(j).cluster+NUM_NODES+1),Inf);
                    end          
                elseif((distance_nodes>threshold)&&(nodes(j).state==1)&&(i~=j))
                    if(distance_nodes<nodes(j).distance)
                        G=addedge(G,i,j,distance_nodes);
                    else
                        G=addedge(G,j,(nodes(j).cluster+NUM_NODES+1),nodes(j).distance);
                    end
                end            
            end
        end
     end
%        for i=1:1:NUM_NODES
%            nodes(i).route = shortestpath(G,i,nodes(i).cluster+NUM_NODES+1);
%        end
end


%      G=graph();
%      for i=1:1:NUM_NODES
%          for j=1:1:NUM_NODES
%              if(source_nodes(i,j)~=0)
%              G=graph(s,t);
%              end
%          end
%          
%      end
% end


% function[go,trace] = distance_based_routing(nodes,matrix,NUM_NODES)
%     global do;
%     global k;
%  i=1;
%  x=1;
%  mat1=triu(matrix);
%  for i=1:NUM_NODES+k
%      for j=1:NUM_NODES+k
%          if(mat1(i,j)~=0)
%              mat(i,j)=x;
%              mat(j,i)=mat(i,j);
%              x=x+1;
%          end
%      end
%  end
% disp(triu(matrix));
% for from=1:NUM_NODES+k        % fill initial matrix with via = to
%     for via=1:NUM_NODES+k
%         for to=1:NUM_NODES+k
%             if(from~=via&&from~=to)
%                 if(via==to&&matrix(from,to)~=0)
%                     go(to,via,from)=matrix(from,to);
%                 else
%                     go(to,via,from)=100;
%                 end
%             else
%                     go(to,via,from)=inf;
%             end
%         end
%     end
% end
% i=0;
% while(i<2)
% for from=1:NUM_NODES+k
%     for to=1:NUM_NODES+k
%         if(from~=to)
%             if(matrix(from,to)~=0)%calculate neighbour node
%                 for x=1:NUM_NODES+k
%                     for y=1:NUM_NODES+k
%                         temp(x,y)=matrix(from,to)+min(go(y,:,to));
%                         if(temp(x,y)<go(y,to,from)&&go(y,to,from)<inf)
%                             go(y,to,from)=temp(x,y);
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% i=i+1;
% end
% disp(go);
% % Specify range
% range = 18;
% trace = zeros;
% for i=1:1:NUM_NODES+k
%     if(nodes(i).distance>range)
%             source = i;
%             dest = nodes(source).cluster+NUM_NODES;
%             trace(1)=source;
%             j=2;
%             while(source~=dest)
%             [row,col]=find(go(dest,:,source)==min(go(dest,:,source)));
%             %q=find(go(dest,:,source)==min(go(dest,:,source)));
%             if (numel(col)>1)
%                 trace(j)=col(1);
%                 source=col(1);
%                 j=j+1;
%             else
%                 trace(j)=col;
%                 source=col;
%                 j=j+1;
%             end
%             end
%             k=1:j-1;
%             disp(trace(k));
%     end
% end
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%