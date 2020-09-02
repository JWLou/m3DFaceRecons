function [f, u0, v0] = updateCamInParam(pts2d, pts3d, R, T)
% Estimate camera intrinsic parameters with known extrinsic parameters and
% the perspective-n-points constraint

nPts = size(pts2d, 1);
pts3d_c = R*pts3d'+repmat(T, [1, nPts]); pts3d_c = pts3d_c'; % camera coordinates
yu = pts2d(:, 1); xu = pts3d_c(:, 1)./pts3d_c(:, 3); % yu(u) = fu*xu(Xc/Zc)+u0
yv = pts2d(:, 2); xv = pts3d_c(:, 2)./pts3d_c(:, 3); % yv(v) = fv*xv(Yc/Zc)+v0

Au = [xu, ones(nPts, 1)]; qu = Au\yu;
Av = [xv, ones(nPts, 1)]; qv = Av\yv;
f = (qu(1)+qv(1))/2;
u0 = qu(2);
v0 = qv(2);

end