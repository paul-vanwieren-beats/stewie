%Gyro orientation demo

%%%
%Comment this line out if you've recalibrated. This calibration was
%performed on headset "2"
%offsets = [1.3867; 1.2012; 1.6681];
%%%

clearvars -except offsets portName
close all;
clc;

%Settings
baudRate = 230400;
timeout = 0.1;
delay = 0.001;
lineWidth = 3;

%Turn off warnings
warning('off');

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

%Set up plot
fig = figure;
% axs = axes;
% hold(axs, 'on');
% grid(axs, 'on');
% axis square;
% view(-62, 10);
% set(axs, 'xlim', [-1.25, 1.25], 'ylim', [-1.25, 1.25], 'zlim', [-1.25, 1.25]);


% fig2 = figure;
% a2 = axes;
% hold(a2, 'on');
% grid(a2, 'on');


%Plot unit circles
% rads = 0 : 0.01 : 2 * pi;
% zero = zeros(1, length(rads));
% cosines = cos(rads);
% sines = sin(rads);
% plot3(axs, cosines, sines, zero, 'r-', 'linewidth', lineWidth);
% plot3(axs, cosines, zero, sines, 'r-', 'linewidth', lineWidth);
% plot3(axs, zero, cosines, sines, 'r-', 'linewidth', lineWidth);

%Pause for figure update
pause(delay);

%Main loop
try
    tStart = now;
    azimuth = 0;
    elevation = 0;
    i = 1;
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
                
                %Get gyroscope values
                gx = hexadecimalToDecimal(fields{6}) / 100000;
                gy = hexadecimalToDecimal(fields{7}) / 100000;
                gz = hexadecimalToDecimal(fields{8}) / 100000;
                
                %Apply offsets
                gx = gx - offsets(1);
                gy = gy - offsets(2);
                gz = gz - offsets(3);
                
                %Integrate
                if exist('timestamp')
                    dt = (now - timestamp) * 24 * 60 * 60;
                    azimuth = azimuth + (gz * dt);
                    elevation = elevation + (gy * dt);
                end
                timestamp = now;
                
                %Delete previous plot
                if exist('p')
                    delete(p);
                end
                
                %Output
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
                [~, ~] = system(sprintf('./sendOSC -h localhost 3456 "%i"', -round(azimuth)));
                
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

%Close port
fprintf(port, '%s\r\n', 'X');
flushinput(port);
fclose(port);

