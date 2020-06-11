classdef Sensor_Node
    %S Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x;
        y;
        G=0;
        type='N';
        battery=0.5;
        sink_x;
        sink_y;
        state = 1;
        cluster;
    end
    
     methods
         function Node = Add_Positions(Node)
%             %S Construct an instance of this class
%             %   Detailed explanation goes here
                 Node.x=rand(1,1)*100;
                 Node.y=rand(1,1)*100;
         end
     end
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
%     end
end

