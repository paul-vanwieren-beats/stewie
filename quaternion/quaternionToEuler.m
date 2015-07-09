function angles = quaternionToEuler(q, eulerPrevDeg)
% Convert quaternion to Euler angles [deg]
%
% Quaternion, q, represents the rotation of the coordinate frame
% from the inertial coordinate frame (Earth) to the local coordinate frame
% (device).
%
% Returned matrix is 3 x N with column k holding the set of Euler angles
% for the kth quaternion, if q is a 4 x N matrix. The 1st row is
% theta (pitch), 2nd row is phi (roll), 3rd row is psi (yaw).
%
% Uses Apple's Euler angle convention: ZXY
%   p_local = R * p_earth    where
%         R = Ry(roll) * Rx(pitch) * Rz(yaw)
%

if (nargin < 2)
	eulerPrevDeg = [0;0;0];
end

if (size(q,1) ~= 4)
	% Transpose so that q is 4 x N matrix where N is the number of
	% quaternions to be converted.
	q = q';
end
N = size(q,2);

angles = zeros(3, N);
for n=1:N
	R = quaternionToRotationMatrix(q(:,n));
	angles(:,n) = [ ...
		asin(R(2,3)); ...
		-atan2(R(1,3), R(3,3)); ...
		-atan2(R(2,1), R(2,2)); ];
end

%% Catch singularity at pitch = +/- 90 deg

% At singularity we can only determine either the sum of
% difference of theta and psi. We use the assumption of
% time continuity to unwrap the two.
% 
% singularity = R(2,3) > 0.99;
% idx = find(singularity);
% if ~isempty(idx)
% 	for n=1:length(idx)
% 		R = quaternionToRotationMatrix(q(:,n));
% 		phiPlusPsi = -atan2(R(3,1), R(3,2));
% 		if (n>1)
% 			phi = angles(2,n-1);
% 		else
% 			phi = eulerPrevDeg(2) * pi/180;
% 		end
% 		angles(:,n) = [ ...
% 			asin(R(2,3));
% 			phi;
% 			mod(phiPlusPsi - phi + pi,2*pi)-pi]; % Wrap to +/- pi
% 	end
% end
% 
% singularity = R(2,3) < -0.99;
% idx = find(singularity);
% if ~isempty(idx)
% 	for n=1:length(idx)
% 		R = quaternionToRotationMatrix(q(:,n));
% 		phiMinusPsi = -atan2(R(3,1), R(3,2));
% 		if (n>1)
% 			phi = angles(2,n-1);
% 		else
% 			phi = eulerPrevDeg(2) * pi/180;
% 		end
% 		angles(:,n) = [ ...
% 			asin(R(2,3));
% 			phi;
% 			mod(phi - phiMinusPsi + pi,2*pi)-pi]; % Wrap to +/- pi
% 	end
% end

% singularity = R(2,3) > 0.99;
% idx = find(singularity);
% if ~isempty(idx)
% 	for n=1:length(idx)
% 		R = quaternionToRotationMatrix(q(:,n));
% 		phiPlusPsi = -atan2(R(3,1), R(3,2));
% 		psi = 0;
% 		phi = phiPlusPsi;
% 		angles(:,n) = [ ...
% 			asin(R(2,3));
% 			phi;
% 			psi];
% 	end
% end
% 
% singularity = R(2,3) < -0.99;
% idx = find(singularity);
% if ~isempty(idx)
% 	for n=1:length(idx)
% 		R = quaternionToRotationMatrix(q(:,n));
% 		phiMinusPsi = -atan2(R(3,1), R(3,2));
% 		phi = 0;
% 		psi = phi-phiMinusPsi;
% 		angles(:,n) = [ ...
% 			asin(R(2,3));
% 			psi;
% 			phi];
% 	end
% end

% Convert to degrees
angles = real(angles) * 180/pi;
