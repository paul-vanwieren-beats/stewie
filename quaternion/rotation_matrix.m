function R = rotation_matrix(axis,angle)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    sa = sin(angle);
    ca = cos(angle);
    switch (axis)
        case 1
            R = [1 0 0; 0 ca sa; 0 -sa ca];
        case 2
            R = [ca 0 -sa; 0 1 0; sa 0 ca];
        case 3
            R = [ca sa 0; -sa ca 0; 0 0 1];
        otherwise
            error('invalid axis %d', axis);
    end
end

