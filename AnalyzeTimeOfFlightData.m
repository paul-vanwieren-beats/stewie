%Analyze time of flight data from sensor array

%Clean up
clear all;
close all;
clc;

%Settings
filePath = 'ToF_Array_Planar_Surface.csv';

%Read data
fid = fopen(filePath);
for i = 1:135
    line = fgetl(fid);
    fields = strsplit(line, ',');
    for j = 2:length(fields)
        raw(i, j - 1) = -str2double(fields{j});
        if raw(i, j - 1) == 0
            raw(i, j - 1) = nan;
        end
    end
    means(i, 1) = mean(raw(i, :));
    stdDevs(i, 1) = std(raw(i, :));
end

%Re-package data into 3d scatter values
x = 49;
y = 44;
for i = 1:135
    xv(i, 1) = x;
    yv(i, 1) = y;
    zv(i, 1) = means(i);
    if x == 0
        x = 49;
        y = y - 5.5;
    else
        x = x - 3.5;
    end
end

%Interpolate data to 1mm spacing
xi = (0:1:49)';
yi = (0:1:44)';
[xg, yg] = meshgrid(xi, yi);
f = scatteredInterpolant(xv, yv, zv);
zg = f(xg, yg);

%Plot data
figure;
axes;
hold on;
grid on;
surf(xg, yg, zg);
plot3(xv, yv, zv, '.', 'markersize', 15);