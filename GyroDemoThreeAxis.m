%3-axis gyro fusion algorithm

%Clean up
clearvars -except offsets portName
close all;
clc;

%Dependencies
addpath(fullfile(cd, 'quaternion'));

%Settings
baudRate = 230400;
timeout = 0.1;
delay = 0.001;

%Gyro biases estimated from calibration procedure
bias = offsets * pi / 180;

%Load gyro data
% t = (0:0.005:10)';
% u = rand(1,3);
% u = u/norm(u);
% raw = bsxfun(@plus,bsxfun(@times,5*sin(2*pi*(t./10)),u),bias); % [dps]

%Algorithm setup
sensorfusion = SensorFusionOIS();
%sensorfusion.bias = (bias' + 0.15*randn(3,1))*degrees;
sensorfusion.bias = bias';

%Open serial port
try
    port = serial(portName);
    set(port, 'BaudRate', baudRate, 'Timeout', timeout);
    fopen(port);
    pause(0.1);
    flushinput(port);
    fprintf(port, '%s\r\n', 'R');
catch err
    error(err.message);
end

%Set up figure
fig = figure;
pause(delay);

%Main loop
try
    tStart = now;
    while true

        %Check for quit condition
        if get(fig, 'currentCharacter') == 'q'
            
            break;
        
        else

            %Get line and parse
            try
                line = fgetl(port);
            catch
                fprintf('fgetl error\n');
                flushinput(port);
                continue;
            end
            fields = splitString(line);
            if strcmp(fields{1}, '4') && (length(fields) == 14)
                
                %Get accelerometer values
                ax = -hexadecimalToDecimal(fields{3}) / 100;
                ay = -hexadecimalToDecimal(fields{4}) / 100;
                az = -hexadecimalToDecimal(fields{5}) / 100;
                
                %Get gyroscope values
                gx = hexadecimalToDecimal(fields{6}) / 100000;
                gy = hexadecimalToDecimal(fields{7}) / 100000;
                gz = hexadecimalToDecimal(fields{8}) / 100000;
                
                %Run algorithm
                t = (now - tStart) * 24 * 60 * 60;
                aRaw = [ax, ay, az];
                gRaw = [gx, gy, gz] * pi / 180;
                sensorfusion.feedAcceleration(aRaw);
                sensorfusion.feedGyro(gRaw, t);
                attitude = sensorfusion.attitude();
                gravity = sensorfusion.gravity()';
                reference = cross([0; 1; 0], gravity);
                reference = reference / norm(reference);
                north = (attitude.dcm())' * reference;
%                 north = north - dot(north, gravity) * gravity;
                heading = atan2d(-north(2), -north(1));
                
%                 fprintf('%0.2f\n', heading);

%                 fprintf('%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n', gravity(1), gravity(2), gravity(3), north(1), north(2), north(3));
                
                %Outputs
                %                 azimuthRadians = azimuth * pi / 180;
                %                 elevationRadians = elevation * pi / 180;
                %                 [x, y, z] = sph2cart(azimuthRadians, elevationRadians, 1);
                %                 p = plot3(axs, [0, x], [0, y], [0, z], 'b-', 'linewidth', lineWidth);
                %
                %                 plot(a2, i, gx, 'k.');
                %                 plot(a2, i, gy, 'b.');
                %                 plot(a2, i, gz, 'r.');
                
                %                 tElapsed = (timestamp - tStart) * 24 * 60 * 60;
                %                 if exist('dt')
                %                     fprintf('%d\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n', i, tElapsed, gx, gy, gz, azimuth, elevation);
                %                 end
                [~, ~] = system(sprintf('./sendOSC -h localhost 3456 "%i"', -round(heading)));
                
                i = i + 1;
                
                flushinput(port);
                
            end
        end
        
        %Pause for keypress
        pause(delay);
        
    end
catch err
    fclose(port);
    for l = 1:length(err.stack)
        disp(['Error in ' err.stack(l).file ' at line ' num2str(err.stack(l).line)]);
    end
    error('Failed');
end

%Close port and close figure
fprintf(port, '%s\r\n', 'X');
flushinput(port);
fclose(port);
close(fig);

%% plots
% figure('Name', mfilename);
% subplot(2,1,1);
% plot(t,w_hat*radians);
% hold on;
% plot(t,raw,'x');
% hold off;
% xlabel('time [sec]');
% ylabel('angular rate [dps]');
% grid on;
% 
% subplot(2,1,2);
% plot(t, q_hat(:,1:3)*2*radians);
% xlabel('time [sec]');
% ylabel('angles [deg]');
% grid on;