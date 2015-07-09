%Calibrate gyros
%Make sure headphones are sitting still on table before starting

clear all;
close all;
clc;

%Settings
portName = '/dev/tty.usbserial-A7042TQ4';
baudRate = 230400;
timeout = 0.1;
delay = 0.001;
n = 100;

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

%Main loop
gxv = nan(n, 1);
gyv = nan(n, 1);
gzv = nan(n, 1);
try
    
    i = 1;
    while i <= n
        
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
            
            %Add values to vectors
            gxv(i, 1) = gx;
            gyv(i, 1) = gy;
            gzv(i, 1) = gz;
            
            fprintf('%0.5f\t%0.5f\t%0.5f\n', gx, gy, gz);
            flushinput(port);
            i = i + 1;
            
        end
    end
    
    %Pause for keypress
    pause(delay);
    
    %Calculate offsets
    offsets = [mean(gxv); mean(gyv); mean(gzv)]
    
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
        
