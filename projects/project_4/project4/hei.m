close all;
clear all;
clc;


for i=1:100
    data = randi(10,10);

    %data
        
    image(data(:,:,1))
    drawnow()
end