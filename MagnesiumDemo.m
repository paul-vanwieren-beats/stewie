%Magnesium orientation demo

clearvars -except offsets sensitivities portName
close all;
clc;

%Settings
baudRate = 230400;
timeout = 0.1;
delay = 0.001;
nPlot = 1;
% lineWidth = 3;

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
% axis square;
% set(axs, 'xlim', [-1.25, 1.25], 'ylim', [-1.25, 1.25]);

%Plot unit circle
% rads = 0:0.01:2*pi;
% unitX = cos(rads);
% unitY = sin(rads);
% plot(axs, unitX, unitY, 'r-', 'linewidth', lineWidth);

%Pause for figure update
% pause(delay);

%Main loop
try
    stop = false;
    i = 1;
    while ~stop
        
        %Check for quit condition
        if get(fig, 'currentCharacter') == 'q'
            stop = true;
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
                ax = hexadecimalToDecimal(fields{3})/100;
                ay = hexadecimalToDecimal(fields{4})/100;
                az = hexadecimalToDecimal(fields{5})/100;
                
                %Get magnetometer values
                mx = hexadecimalToDecimal(fields{12})/100;
                my = hexadecimalToDecimal(fields{13})/100;
                mz = hexadecimalToDecimal(fields{14})/100;
                
                %Create sensor vectors
                a = [ax; ay; az];
                m = [mx; my; mz];
                
                %Apply magnetometer offsets
                m = (m - offsets);
                
                %Create common axes
                a = [a(2); a(1); a(3)];
                m = [-m(2); -m(1); -m(3)];
                
                %Convert to x/y angle
                thetaRad = atan2(m(2), m(1));
                thetaDeg = thetaRad*180/pi;
                
                %Convert to x and y on unit circle
                x = cos(thetaRad);
                y = sin(thetaRad);
                
                %Check for plot condition
                if i == nPlot
                    
                    %Delete previous plot
                    %                     if exist('p')
                    %                         delete(p);
                    %                     end
                    
                    %Output
                    %                     p = plot(axs, [0, x], [0, y], 'b-', 'linewidth', lineWidth);
%                                          fprintf('%0.2f\t%0.2f\t%0.2f\t%0.2f\n', m(1), m(2), m(3), thetaDeg);
                    [~, ~] = system(sprintf('./sendOSC -h localhost 3456 "%i"', -round(thetaDeg)));
                    flushinput(port);
                    
                    %Reset
                    i = 1;
                    
                else
                    i = i + 1;
                end
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
        
