%Magnesium orientation demo

clear all;
close all;
clc;

%Settings
portName = '/dev/tty.usbserial-A7042TQI';
baudRate = 230400;
timeout = 0.1;
delay = 0.001;
nPlot = 20;
nArray = 600; %60 seconds of data
aLimit = 75;

%Turn off warnings
warning('off');

%Create arrays
x = nan(nArray, 1);
y = nan(nArray, 1);
z = nan(nArray, 1);

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
axs = axes;
hold(axs, 'on');
axis square;
view(45, 15);
% xlim([-aLimit, aLimit]);
% ylim([-aLimit, aLimit]);
% zlim([-aLimit, aLimit]);
grid on;
set(axs, 'fontsize', 20);
xlabel('x', 'fontsize', 20);
ylabel('y', 'fontsize', 20);
zlabel('z', 'fontsize', 20);

%Pause for figure update
pause(delay);

%Main loop
try
    stop = false;
    i = 1;
    n = 0;
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
                
                %Get magnetometer values
                mx = hexadecimalToDecimal(fields{12})/100;
                my = hexadecimalToDecimal(fields{13})/100;
                mz = hexadecimalToDecimal(fields{14})/100;

                %Check for plot condition
                if i == nPlot
                    
                    %Update arrays
                    x = circshift(x, 1);
                    x(1) = mx;
                    y = circshift(y, 1);
                    y(1) = my;
                    z = circshift(z, 1);
                    z(1) = mz;
                    n = n + 1;
                    
                    %Calculate offsets and sensitivities
                    if n > length(x)
                        xc = x(1:end);
                        yc = y(1:end);
                        zc = z(1:end);
                    else
                        xc = x(1:n);
                        yc = y(1:n);
                        zc = z(1:n);
                    end
                    xOffset = mean([min(xc), max(xc)]);
                    yOffset = mean([min(yc), max(yc)]);
                    zOffset = mean([min(zc), max(zc)]);
                    xSensitivity = (max(xc) - min(xc))/2;
                    ySensitivity = (max(yc) - min(yc))/2;
                    zSensitivity = (max(zc) - min(zc))/2;
                    offsets = [xOffset; yOffset; zOffset];
                    sensitivities = [xSensitivity; ySensitivity; zSensitivity];
                    
                    %Output
                    fprintf('%d\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n', n, offsets(1), offsets(2), offsets(3), sensitivities(1), sensitivities(2), sensitivities(3));
                    scatter3(axs, mx, my, mz, 'ro', 'markerfacecolor', 'r');
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
        
