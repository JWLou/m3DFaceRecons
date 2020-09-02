function params = estCamParam(params, lan2d)
% Estimate initial camera parameters according to the landmark constraint

% Detected 2D landmarks on the image
nLan = size(lan2d, 1);
lan2d_h = [lan2d, ones(nLan, 1)];
mu_shape3d = params.mu_id + params.mu_exp;
mu_shape3d = reshape(mu_shape3d, [3, length(mu_shape3d)/3]);
lan3d = mu_shape3d(:, params.idxKpts); lan3d = lan3d';

% Estimate camera parameters using EPnP 
A = [params.f, 0, params.u0; 0, params.f, params.v0; 0, 0, 1]; % the initial camera intrinsic matrix 
lan3d_h = [lan3d, ones(nLan, 1)];
[R, T, ~, ~] = EPnP(lan2d_h, lan3d_h, A); % estimate camera extrinsic parameters
[xAngle, yAngle, zAngle] = dcm2angle(R, 'XYZ'); % get rotation angles from the rotation matrix, XYZ-extrinsic, left-handed
params.vecAngle = [-xAngle; -yAngle; -zAngle]; % extrinsic ;
params.T = double(T);

end
