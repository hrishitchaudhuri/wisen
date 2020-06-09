clc; clear; close;

%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%

x_max = 100;
y_max = 100;

sink.x = 0.5 * x_max;
sink.y = 0.5 * y_max;

NUM_NODES = 100;

Eo=0.5;

ETX=50*0.0000001;
ERX=50*0.0000001;

Efs=10*0.0000000001;
Emp=0.0013*0.0000000001;

rounds = 100;

%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%

do=sqrt(Efs/Emp);
dead_nodes = 0;

figure(1);

for i = 1:1:NUM_NODES
    S(i).x = rand(1, 1) * x_max;
    S(i).y = rand(1, 1) * y_max;
    S(i).battery = Eo;
    S(i).state = 1;
    plot(S(i).x, S(i).y, 'black o');
    hold on;
end

plot(sink.x, sink.y, 'green x', "MarkerSize", 12);
hold on;

for i = 1:1:rounds
    for i = 1:1:NUM_NODES
        if (S(i).state == 1)
            
            distance=sqrt((S(i).x - sink.x)^2 + (S(i).y - sink.y)^2);
            if (distance>do)
                S(i).battery = S(i).battery - ((ETX)*(4096) + Emp*4096*(distance^4)); 
            end
        
            if (distance<=do)
                S(i).battery = S(i).battery - ((ETX)*(4096) + Efs*4096*(distance^2)); 
            end
        
            if (S(i).battery <= 0)
                dead_nodes = dead_nodes + 1;
                S(i).state = 0;
                plot(S(i).x, S(i).y, 'black x');
            end
            
        end
    end
    
    pause(0.5);
end