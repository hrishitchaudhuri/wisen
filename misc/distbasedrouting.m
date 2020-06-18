clc;
clear all;
close all;
n=6;
over=0;
for i=1:n    %initialize
    for j=i:n       
        if(i==j)
            matrix(i,j)=0;
        else
            matrix(i,j)=randi(9,1,1);
            matrix(j,i)=matrix(i,j);
        end
    end
end
disp(matrix);
disp("matrix");
 i=1;
 x=1;
 mat1=triu(matrix);
 %disp(mat1);
 for i=1:n
     for j=1:n
         if(mat1(i,j)~=0)
             mat(i,j)=x;
             mat(j,i)=mat(i,j);
             x=x+1;      
         end
     end
 end
disp(mat1);
disp('mat1');
disp(mat);
disp('mat');

for from=1:n    % fill initial matrix with via = to
    for via=1:n
        for to=1:n           
            if(from~=via&&from~=to)
                if(via==to&&matrix(from,to)~=0)
                    go(to,via,from)=matrix(from,to);
                else
                    go(to,via,from)=100;
                end
            else
                    go(to,via,from)=inf;
            end           
        end
    end
end
disp(go);
disp('go');

i=0;
while(i<2)
for from=1:n
    for to=1:n
        if(from~=to)
            if(matrix(from,to)~=0)%calculate neighbour node
                for x=1:n
                    for y=1:n                    
                        temp(x,y)=matrix(from,to)+min(go(y,:,to));
                        if(temp(x,y)<go(y,to,from)&&go(y,to,from)<inf)
                            go(y,to,from)=temp(x,y);
                        end                      
                    end
                end
            end
        end
    end
end
i=i+1;
end
disp(go);
disp('go final');
source=input('Enter the source node: ');
dest=input('Enter the destination node: ');
trace(1)=source;
j=2;
% [row,col]=find(go(dest,:,source)==min(go(dest,:,source)));
while(source~=dest)
[row,col]=find(go(dest,:,source)==min(go(dest,:,source)));
%q=find(go(dest,:,source)==min(go(dest,:,source)));
if (numel(col)>1)
    trace(j)=col(1);
    source=col(1);
    j=j+1;
else
    trace(j)=col;
    source=col;
    j=j+1;
end
end
k=1:j-1;
disp(trace(k));
