classdef Sensor_Node
    %S Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xd;
        yd;
        G=0;
        type='N';
        E=0.5;
        base_dis;
        cluster;
    end
    
     methods
         function Node = Add_Positions(Node)
%             %S Construct an instance of this class
%             %   Detailed explanation goes here
                 Node.xd=rand(1,1)*100;
                 Node.yd=rand(1,1)*100;
         end
     end
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
%     end
end

