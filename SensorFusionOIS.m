classdef SensorFusionOIS < handle
    % Three axis sensor fusion (gyro rate integration only)
    % Adam Howell (adam_howell@apple.com)
    
    properties
        q = Quaternion();
        time = [];
        delta = 0;
        
        wc = zeros(3,1);
        bias = [];
        
        filteredGyro = [];
        gyroLpf = exp(-2*pi*30/200);
        
        filteredAccel = [];
        accelLpf = exp(-2*pi*1/100);
    end
    
    properties(Dependent)
       attitude;
       rate;
       gravity;
    end
    
    methods
        function SF = SensorFusionOIS()
            
        end
        
        function qo = get.attitude(SF)
            qo = SF.q;
        end
        
        function wo = get.rate(SF)
            wo = SF.wc;
        end
        
        function go = get.gravity(SF)
            go = SF.filteredAccel/norm(SF.filteredAccel);
        end
        
        function feedGyro(SF,w,timestamp)
            
            % assume that bias is available prior to 
            if isempty(SF.bias), return; end
            
            % correct for bias
            SF.wc = w - SF.bias;
            
            % LPF for prediction
            if isempty(SF.filteredGyro)
                SF.filteredGyro = SF.wc;
            else
                SF.filteredGyro = SF.gyroLpf*SF.filteredGyro + (1-SF.gyroLpf)*SF.wc;
            end
            
            % integrate the state
            SF.q = SF.predict(timestamp, SF.wc);
            if isempty(SF.time)
                SF.delta = 0;
            else
                SF.delta = timestamp - SF.time;
            end
            SF.time = timestamp;
        end
        
        function feedAcceleration(SF,a,~)
            % LPF for gravity unit vector
            if isempty(SF.filteredAccel)
                SF.filteredAccel = a;
            else
                SF.filteredAccel = SF.accelLpf*SF.filteredAccel + (1-SF.accelLpf)*a;
            end
        end
        
        function feedGyroBias(SF, ~, event)
            % event handler
            SF.bias = event.bias;
        end
        
        function qp = predict(SF,timestamp,w)
            
            if ~exist('w','var'), w = SF.filteredGyro; end
            
            if isempty(SF.time) || (SF.time == timestamp),
                qp = SF.q;
            else
                % integrate the state
                dt = timestamp - SF.time;
                qp = Quaternion(w*dt)*SF.q;
            end
            
        end
        
    end
    
end

