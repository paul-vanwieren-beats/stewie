%3-axis gyro fusion algorithm

%Dependencies
addpath(fullfile(cd, 'quaternion'));

%% constants
degrees = pi/180;
radians = 1/degrees;

%% gyro bias estimated from registration period
bias = 3.5*randn(1,3);

%% load gyro data
t = (0:0.005:10)';
u = rand(1,3);
u = u/norm(u);
raw = bsxfun(@plus,bsxfun(@times,5*sin(2*pi*(t./10)),u),bias); % [dps]

%% algorithm setup
sensorfusion = SensorFusionOIS();
sensorfusion.bias = (bias' + 0.15*randn(3,1))*degrees;

%% loop it baby
N = numel(t);
w_hat = zeros(N,3);
q_hat = zeros(N,4);
for k = 1:N,   
    w = raw(k,:)'*pi/180; % [rad/sec]
    
    sensorfusion.feedGyro(w, t(k));
    
    q_hat(k,:) = sensorfusion.attitude()';
    w_hat(k,:) = sensorfusion.rate()';
end

%% plots
figure('Name', mfilename);
subplot(2,1,1);
plot(t,w_hat*radians);
hold on;
plot(t,raw,'x');
hold off;
xlabel('time [sec]');
ylabel('angular rate [dps]');
grid on;

subplot(2,1,2);
plot(t, q_hat(:,1:3)*2*radians);
xlabel('time [sec]');
ylabel('angles [deg]');
grid on;