%Test 3d plot

%Clean up
clear all;
close all;
clc;

n = 1000;
delay = 0.001;

figure;
axes;
hold on;
x = nan(n, 1);
y = nan(n, 1);
z = nan(n, 1);
for i = 1:n
    
    %Random numbers
    xi = 1 - 2*rand;
    yi = 1 - 2*rand;
    zi = 1 - 2*rand;
    
    %Scale to unit sphere
    scale = 1/sqrt(xi^2 + yi^2 + zi^2);
    xi = xi*scale;
    yi = yi*scale;
    zi = zi*scale;
    
    %Save in array
    x(i) = xi;
    y(i) = yi;
    z(i) = zi;
    
end

%Plot
scatter3(x, y, z, 'ro');
axis square;