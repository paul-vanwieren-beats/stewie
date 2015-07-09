function q = quaternionFromRotationMatrix(m)

q = zeros(1, 4);


if trace(m) > 0.0
    t = m(1, 1) + m(2, 2) + m(3, 3) + 1.0;
    s = 0.5 / sqrt(t);
    q(4) = s * t;
    q(3) = ( m(1, 2) - m(2, 1) ) * s; 
    q(2) = ( m(3, 1) - m(1, 3) ) * s; 
    q(1) = ( m(2, 3) - m(3, 2) ) * s;
elseif m(1, 1) > m(2, 2) && m(1, 1) > m(3, 3)
    t = m(1, 1) - m(2, 2) - m(3, 3) + 1.0;
    s = 0.5 / sqrt(t);
    q(1) = s * t;
    q(2) = ( m(1, 2) + m(2, 1) ) * s;
    q(3) = ( m(3, 1) + m(1, 3) ) * s; 
    q(4) = ( m(2, 3) - m(3, 2) ) * s;
elseif m(2, 2) > m(3, 3)
    t = -m(1, 1) + m(2, 2) - m(3, 3) + 1.0;
    s = 0.5 / sqrt(t);
    q(2) = s * t; 
    q(1) = ( m(1, 2) + m(2, 1) ) * s;
    q(4) = ( m(3, 1) - m(1, 3) ) * s;
    q(3) = ( m(2, 3) + m(3, 2) ) * s;
else
    t = -m(1, 1) - m(2, 2) + m(3, 3) + 1.0;
    s = 0.5 / sqrt(t);
    q(3) = s * t;
    q(4) = ( m(1, 2) - m(2, 1) ) * s;
    q(1) = ( m(3, 1) + m(1, 3) ) * s;
    q(2) = ( m(2, 3) + m(3, 2) ) * s;
end

% Calculate the trace of the matrix T from the equation:
% 
%                 2     2     2
%       T = 4 - 4x  - 4y  - 4z
% 
%                  2    2    2
%         = 4( 1 -x  - y  - z )
% 
%         = 1 + mat(0) + mat(5) + mat(10)
% 
% 
%     If the trace of the matrix is greater than zero, then
%     perform an "instant" calculation.
%     Important note wrt. rouning errors:
% 
%     Test if ( T > 0.00000001 ) to avoid large distortions!
% 
%       S = sqrt(T) * 2;
%       X = ( mat(9) - mat(6) ) / S;
%       Y = ( mat(2) - mat(8) ) / S;
%       Z = ( mat(4) - mat(1) ) / S;
%       W = 0.25 * S;
% 
%     If the trace of the matrix is equal to zero then identify
%     which major diagonal element has the greatest value.
%     Depending on this, calculate the following:
% 
%     if ( mat(0) > mat(5) && mat(0) > mat(10) )  {	// Column 0: 
%         S  = sqrt( 1.0 + mat(0) - mat(5) - mat(10) ) * 2;
%         X = 0.25 * S;
%         Y = (mat(4) + mat(1) ) / S;
%         Z = (mat(2) + mat(8) ) / S;
%         W = (mat(9) - mat(6) ) / S;
%     } else if ( mat(5) > mat(10) ) {			// Column 1: 
%         S  = sqrt( 1.0 + mat(5) - mat(0) - mat(10) ) * 2;
%         X = (mat(4) + mat(1) ) / S;
%         Y = 0.25 * S;
%         Z = (mat(9) + mat(6) ) / S;
%         W = (mat(2) - mat(8) ) / S;
%     } else {						// Column 2:
%         S  = sqrt( 1.0 + mat(10) - mat(0) - mat(5) ) * 2;
%         X = (mat(2) + mat(8) ) / S;
%         Y = (mat(9) + mat(6) ) / S;
%         Z = 0.25 * S;
%         W = (mat(4) - mat(1) ) / S;
%     }
% 
