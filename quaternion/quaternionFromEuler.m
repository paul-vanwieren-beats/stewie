function qEuler = quaternionFromEuler(pitch,roll,yaw)
% Convert Euler angles [deg] to quaternion representation.
%
% Quaternion, q, represents the rotation of the coordinate frame
% from the inertial coordinate frame (Earth) to the local coordinate frame
% (device).
%
% Returned matrix is 4 x N with column k holding the quaternion for the
% kth set of Euler angles.
%
% Uses Apple's Euler angle convention: ZXY
%   p_local = R * p_earth    where
%         R = Ry(roll) * Rx(pitch) * Rz(yaw)
%
% Call syntax:
% 1.  quaternionFromEuler(angleMatrix)
%        where angleMatrix is 3xN consisting of N Euler angles.
%        The 1st row is pitch, 2nd row is roll and 3rd row is yaw.
% 2.  quaternionFromEuler(pitch, roll, yaw)
%
%
% Ranges:
%    -180 <= roll < 180

if nargin == 1
	if size(pitch,1) == 3
		a = pitch;
	else
		error('Invalid arguments');
	end
else
	if length(pitch) ~= length(roll) || length(pitch) ~= length(yaw)
		error('input lengths not consistent');
	end
	% Make into row vectors
	if (size(pitch,1) > 1), pitch = pitch'; end
	if (size(roll,1) > 1), roll = roll'; end
	if (size(yaw,1) > 1), yaw = yaw'; end

	a = [pitch;roll;yaw];
end

% Convert from [deg] to [rad] and change to half angles
a = 0.5 * a * pi/180;

cosA = cos(a);
sinA = sin(a);

theta2 = 1;
phi2   = 2;
psi2   = 3;

% > CodeGeneration[Matlab](qEuler(2*theta2, 2*phi2, 2*psi2), resultname = "qEuler");
qEuler = [-sinA(psi2,:) .* cosA(theta2,:) .* sinA(phi2,:) + cosA(psi2,:) .* sinA(theta2,:) .* cosA(phi2,:); cosA(psi2,:) .* cosA(theta2,:) .* sinA(phi2,:) + sinA(psi2,:) .* sinA(theta2,:) .* cosA(phi2,:); cosA(psi2,:) .* sinA(theta2,:) .* sinA(phi2,:) + sinA(psi2,:) .* cosA(theta2,:) .* cosA(phi2,:); -sinA(psi2,:) .* sinA(theta2,:) .* sinA(phi2,:) + cosA(psi2,:) .* cosA(theta2,:) .* cosA(phi2,:);];

% BKBK: Don't allow q(4,:) to be negative, this will keep cos(alpha) unique
% by restricting them to quadrants I and IV.
% idxNeg = qEuler(4,:) < 0;
% qEuler(:,idxNeg) = -qEuler(:,idxNeg);
