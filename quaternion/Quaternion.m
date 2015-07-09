classdef Quaternion
    %Quaternion Quaternion class
    %   This class implements quaternions consistent with the CoreMotion
    %   representation for attitude.
    %
    %   Features:
    % - Constructor can be supplied rotation matrices, euler angles and angle axis forms are implicit
    % - Static methods for conversions between rotation matrices, euler angles and angle axis forms
    % - Conjugation via ~
    % - Quaternion multiplication for quaternions or {Nx4,4xN} matrices
    % - Indexing via .{x,y,z,w}
    % - Rotation of vector by a quaternion via .rotate()
    % - Interpolation using SLERP via Quaternion.slerp()
    
    properties(Access=private)
        q;
    end
    
	properties(Dependent=true)
		gravity;
		north;
	end

    methods
        
        function Q = Quaternion(in, varargin)
            % constructor can be supplied with
            % - no argument => defaults to [0; 0; 0; 1]
            % - quaternion as a class or 4D vector/matrix
            % - DCM
            % - angle-axis 3d vector
            % - euler angles as 3d vector plus an ordering as a string (i.e. '321')
            if ~exist('in', 'var'), in = [0; 0; 0; 1]; end
            
            [m,n] = size(in);
            if ((m*n > 4) && ((m == 4)||(n == 4)))
                % recursive constructor for an array of Quaternions from
                % 4xN or Nx4 matrix
                if ((n == 4)&&(m ~= 4)), in = in'; end
                Q = [];
                for i = 1:m,
                    Q = [Q Quaternion(in(:,i))];
                end
            else
                % construction of a single quaternion
                if isa(in, 'Quaternion')
                    Q.q = in.q;
                elseif isnumeric(in)
                    switch(numel(in))
                        case 4
                            % quaternion as vector
                            Q.q = in;
                        case 9
                            % rotation matrix
							Q = Quaternion.dcm(in);
                        case 3
                            if (nargin == 2) && isa(varargin{1}, 'char')
                                % euler angles
                                Q = Quaternion.euler(in, varargin{1});
                            else
                                % axis-angle
                                Q = Quaternion.aaxis(in);
                            end
                        otherwise
                            error('unknown conversion for input of length %d', numel(in));
                    end
                else
                    error('unsupported type');
                end
            end
        end
        
        function Q = set.q(Q, in)
            if ~isa(in, 'double')
                error('values must be of type double');
            end
            
            % normalize and set constant to be positive
            Q.q = in(:)/norm(in);
            %             if (Q.q(4) < 0), Q.q = -Q.q; end
        end
        
		function n = get.north(Q)
			n = [2*(Q.q(4)*Q.q(4) + Q.q(1)*Q.q(1)) - 1;
				2*(Q.q(1)*Q.q(2) - Q.q(4)*Q.q(3));
				2*(Q.q(1)*Q.q(3) + Q.q(4)*Q.q(2))];
		end

		function Q = set.north(Q,~)
			error('shouldn''t be setting north alone, used northAndGravity');
		end

		function g = get.gravity(Q)
			g = [2*(-Q.q(1)*Q.q(3) + Q.q(4)*Q.q(2));
				2*(-Q.q(2)*Q.q(3) - Q.q(4)*Q.q(1));
				-1 + 2*(Q.q(1)*Q.q(1) + Q.q(2)*Q.q(2))];
		end

		function Q = set.gravity(Q,g)
			len = norm(g);
			Q.q(1) = -(g(2))/len;
			Q.q(2) = (g(1))/len;
			Q.q(3) = 0;
			Q.q(4) = sqrt((1 - (g(3))/len) * 0.5); % from half-angle identity

			lenSquared = Q.q(1)*Q.q(1) + Q.q(2)*Q.q(2);

			if (lenSquared > 1e-06)
				scale = sqrt((1 - Q.q(4)*Q.q(4))/lenSquared);
				Q.q(1) = Q.q(1)*scale;
				Q.q(2) = Q.q(2)*scale;
			else
				q = zeros(4,1);
				% If device is approximately upside down, snap to that orientation
				if (g(3) > 0)
					Q.q(1) = 1;
				else	% else, we're almost exactly right-side up (Z vertical), to snap to identity
					Q.q(4) = 1;
				end
			end

		end

		

        function v = double(Q)
            v = [Q.q];
        end
        
        function disp(Q)
            if isempty(Q), return; end
            disp([Q.q]);
        end
        
        function v = subsref(Q,s)
            switch(s(1).type)
                case '()'
                    v = Q(s.subs{:});
                case '.'
                    elements = {'x','y','z','w'};
                    switch s(1).subs
                        case elements
                            idx = strcmp(s(1).subs,elements);
                            v = Q.q(idx);
                        case 'disp'
                            disp(Q)
                        case {'dcm','euler','aaxis'}
                            v = Quaternion.(s(1).subs)(Q);
                        otherwise
                            if length(s)>1
                                v = Q.(s(1).subs)(s(2).subs{:});
                            else
                                v = Q.(s.subs);
                            end
                    end
                otherwise
            end
        end
        
        function M = lmatrix(Q)
            % matrix form for left multiplication of quaternions
            % ql*qr = M(ql)*qr
            %
            % M = q(4)*I + [-skew(dq) dq; -dq' 0] where dq = q(1:3)
            %
            M = [  Q.q(4)  Q.q(3) -Q.q(2) Q.q(1); ...
                  -Q.q(3)  Q.q(4)  Q.q(1) Q.q(2); ...
                   Q.q(2) -Q.q(1)  Q.q(4) Q.q(3); ...
                  -Q.q(1) -Q.q(2) -Q.q(3) Q.q(4) ];
        end
        
        function M = rmatrix(Q)
            % matrix form for right multiplication of quaternions
            % ql*qr = M(qr)*ql
            %
            % M = q(4)*I + [skew(dq) dq; -dq' 0] where dq = q(1:3)
            %
            M = [  Q.q(4) -Q.q(3)  Q.q(2) Q.q(1); ...
                   Q.q(3)  Q.q(4) -Q.q(1) Q.q(2); ...
                  -Q.q(2)  Q.q(1)  Q.q(4) Q.q(3); ...
                  -Q.q(1) -Q.q(2) -Q.q(3) Q.q(4) ];
        end
        
        function vo = rotate(Q, vi)
            % rotation of the input which can be either a quaternion or vector
            %
            
            if isa(vi, 'Quaternion')
                % equivalent to Q.lmatrix()*(~Q).rmatrix()
                vo = Quaternion([Q.dcm() zeros(3,1); zeros(1,3) 1]*vi.q);
            elseif (numel(vi) == 3)
                vo = Q.dcm()*vi;
            else
                error('incorrect dimension %d', numel(vi));
            end
        end
        
        
        %% overloaded operators
        function r = not(Q)
            % quaternion conjugate
            r = Quaternion([-Q.q(1:3); Q.q(4)]);
        end
        
        function r = mtimes(Ql, Qr)
            % multiply two quaternions such that the order of rotation is
            % consistent with rotation matrices
            % i.e. ql*qr = dcm(ql)*dcm(qr)
            %
            % see rightMultiplyQuaternions.m
            
            if ~isa(Ql, 'Quaternion')
                Ql = Quaternion(Ql);
            end
            
            if ~isa(Qr, 'Quaternion')
                Qr = Quaternion(Qr);
            end
            
            if (numel(Ql) == 1)
                r = Quaternion(Ql.lmatrix()*double(Qr));
            elseif (numel(Qr) == 1)
                r = Quaternion(Qr.rmatrix()*double(Ql));
            else
                error('Quaternion.mtimes undefined for two Quaternion arrays');
            end
        end
        
    end
    
    methods(Static)
        
        function [Qo,P] = northAndGravity(ni,gi,w)
				if ~exist('w','var'), w = ones(1,2); end
				[qi, P] = wahba([ni; gi],[1 0 0; 0 0 -1], w, 'foam2');
                Qo = Quaternion(qi);
        end
        
        function out = dcm(in)
            % conversion to/from DCM
            if isa(in, 'Quaternion')
                % convert from a quaternion to a DCM
                out = quaternionToRotationMatrix(in.q);
            elseif (numel(in) == 9)
                % convert from a DCM to a quaternion
                if (abs(det(in)) >= 1 + 10*eps),
                    warning('dcm is not normal!');
                end
                out = Quaternion(quaternionFromRotationMatrix(in)');
            end
        end
        
        function out = euler(in, order)
            % conversion to/from Euler angles (Apple form)
            
            if exist('order','var') && ~strcmp(order, '321'),
                error('order %s not currently implemented');
            end
            
            if isa(in, 'Quaternion')
                % convert from a quaternion to euler angles
                out = quaternionToEuler(in.q);
            elseif (numel(in) == 3)
                % convert from a DCM to a quaternion
                out = Quaternion(quaternionFromEuler(in));
            end
        end
        
        function out = aaxis(in)
            % conversion to/from axis/angle
            if isa(in, 'Quaternion')
                % convert from a quaternion to a angle-axis
                ndq = norm(in.q(1:3));
                theta = 2*atan2(ndq,in.q(4));
                out = zeros(3,1);
                if (ndq > eps)
                    out = theta*in.q(1:3)/ndq;
                end
            elseif (numel(in) == 3)
                % convert from an angle-axis to a quaternion
                theta = norm(in);
                u = zeros(3,1);
                if (theta > eps)
                    u = in(:)/theta;
                end
                out = Quaternion([sin(0.5*theta)*u; cos(0.5*theta)]);
            else
                error('argument must be a Quaternion or 3x1');
            end
        end
        
        function q = slerp(ti, qi, t)
            % spherical linear interpolation (SLERP) of quaternions
            % this ensures that the interpolated value is also a quaternion
            % and linear interpolates through the angle of rotation
            
            if ~isa(qi, 'Quaternion')
                qi = Quaternion(qi);
            end
            
            for k = 1:numel(t),
                % determine the time indices to interpolate between
                exact = ti == t(k);
                if any(exact),
                    q(k) = qi(exact);
                    continue;
                end
                
                below = find(ti < t(k), 1, 'last');
                if isempty(below)
                    error('query point below range of interpolant (%f < %f)', ...
                        t(k), ti(1));
                elseif (below == numel(ti))
                    error('query point above range of interpolant (%f > %f)', ...
                        t(k), ti(end));
                end
                above = below+1;
                
                % now do SLERP
                theta = (t(k)-ti(below))/(ti(above)-ti(below));
                aax = Quaternion.aaxis(qi(above)*(~qi(below)));
                q(k) = Quaternion(theta*aax)*qi(below);
            end
        end
        
        function test_slerp()
            % unit test for SLERP algorithm
            N = 10;
            step = 2*pi/N;
            axs = rand(3,1);
            dq = Quaternion(step*axs/norm(axs));
            
            qi = Quaternion();
            for n = 1:N,
                qi = [qi dq*qi(end)];
            end
            ti = 0:N;
            
            figure('Name', 'Quaternion: test_slerp');
            plot(ti,qi,'o-');
            
            for t = N*rand(1,20),
                qo = Quaternion.slerp(ti,qi,t);
                hold on;
                plot(t*ones(4,1),double(qo),'kx:');
                hold off;
            end
        end

		function test_initialize
			R = so3_exp(2*pi*rand(3,1));
			q = Quaternion(R);

			noise = [0.2/(3*40) 0.2/3];
			inertial = [1 0 0; 0 0 -1]; % MN and gravity

			N = 1;
			w = [];
			r = [];
			x = [];
			for k = 1:size(inertial,1),
				w = [w; 1/noise(k)*ones(N,1)];
				r = [r; repmat(inertial(k,:), N, 1)];
				x = [x; bsxfun(@plus, noise(k)*randn(N,3), inertial(k,:)*R')];
			end
% 			x(1,:) = [1 0 0];
			q_ng = Quaternion();
			[q_ng, P] = q_ng.northAndGravity(x(1,:),x(2,:),w);
			dq_ng = double(q_ng*(~q));
            P
			q_hat = Quaternion();
			q_hat.gravity = x(2,:);
            dq = double(q_hat*(~q));
            
            fprintf('estimation error:\n');
            fprintf('\tnorthAndGravity:\t'); fprintf('%f ', 2*asind(dq_ng(1:3))); fprintf(' = %f\n', norm(2*asind(dq_ng(1:3))));
            fprintf('\t\t\t\t'); fprintf('%f ', 2*asind(sqrt(diag(P)))); fprintf(' = %f\n', norm(2*asind(sqrt(diag(P)))));
            fprintf('\tgravity:\t\t'); fprintf('%f ', 2*asind(dq(1:3))); fprintf(' = %f\n', norm(2*asind(dq(1:3))));
			
            [double(q)';
            double(q_ng)';
            double(q_hat)']
            
            x
			x_hat = [q_ng.north'; q_ng.gravity']
			x_hat = [q_hat.north'; q_hat.gravity']

		end
    end
end

