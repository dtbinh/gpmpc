clear; clc; close all;
for i=1:200
    y1(i) = windspeed(i,100,120, 'sine');
    y2(i) = 0.1*windspeed(i,80,130, 'step');
    y3(i) = 0.1*windspeed(i,1,200, 'sine');
end

figure;
subplot(3,1,1); plot(y1); ylabel('X');
subplot(3,1,2); plot(y2); ylabel('Y');
subplot(3,1,3); plot(y3); ylabel('Z'); xlabel('Times');
