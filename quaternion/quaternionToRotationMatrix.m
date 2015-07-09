function M = quaternionToRotationMatrix(q)
% Convert quaternion, q, to rotation matrix, M, where M transforms
% points from the a reference coordinate frame to a new coordinate frame
% that is specifed as a rotation by q from the reference coordinate frame.
%
% For our case of interest, we typically define M as the rotation matrix
% to transform the inertial coordinate frame (Earth) to the local
% coordinate frame (device), where q is attitude quaternion for the device.
%
% pLocal = M  pInertial
%        = q* pInertial q
%
% If q consists of N quaternions represented as a 4 x N matrix, then
% M will be 3 x 3 x N.

if (size(q,1) ~= 4)
	% Transpose so that q is 4 x N matrix where N is the number of
	% quaternions to be converted into rotation matrices.
	q = q';
end
N = size(q,2);

% Calculate rotation matrix for each quaternion
x2 = q(1,:) + q(1,:);
y2 = q(2,:) + q(2,:);
z2 = q(3,:) + q(3,:);
xx2 = q(1,:) .* x2;
yy2 = q(2,:) .* y2;
zz2 = q(3,:) .* z2;
yz2 = q(2,:) .* z2;
wx2 = q(4,:) .* x2;
xy2 = q(1,:) .* y2;
wz2 = q(4,:) .* z2;
xz2 = q(1,:) .* z2;
wy2 = q(4,:) .* y2;

M = zeros(3,3,N);
for n=1:N
    M(:,:,n) = [ ...
        1.0 - yy2(n) - zz2(n), xy2(n) + wz2(n),       xz2(n) - wy2(n); ...
        xy2(n) - wz2(n),       1.0 - xx2(n) - zz2(n), yz2(n) + wx2(n); ...
        xz2(n) + wy2(n),       yz2(n) - wx2(n),       1.0 - xx2(n) - yy2(n) ...
    ];
end
    
% 	matrix(1, 1) = 1.0 - yy2 - zz2;
% 	matrix(2, 2) = 1.0 - xx2 - zz2;
% 	matrix(3, 3) = 1.0 - xx2 - yy2;
% 	matrix(3, 2) =	yz2 - wx2;
% 	matrix(2, 3) =	yz2 + wx2;
% 	matrix(2, 1) =	xy2 - wz2;
% 	matrix(1, 2) =	xy2 + wz2;
% 	matrix(1, 3) =	xz2 - wy2;
% 	matrix(3, 1) =	xz2 + wy2;
