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

% Creating a random sensor network
figure(1);

for i = 1:1:NUM_NODES
    S(i, 1) = rand(1, 1) * x_max;
    S(i, 2) = rand(1, 1) * y_max;
    
    plot(S(i, 1), S(i, 2), 'red o');
    hold on;
end

disp(S);

plot(sink.x, sink.y, 'black x');
hold on;

[idx, C] = kmeans(S, k);
disp(idx);
disp(C);
plot(C(:,1), C(:,2), 'black x', 'MarkerSize', 15);
hold on;