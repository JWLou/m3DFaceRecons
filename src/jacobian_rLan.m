function [J_rLan, rLan] = jacobian_rLan(params, lan2d)
% Calculate the Jacobian matrix of landmark residuals with respect to the
% parameters based on IRLS solver - [alpha, beta, delta, vecAngle(phi(angleX), 
% theta(angleY), psi(angleZ)), T, f, (u0, v0), gamma(3*9)]

%% Parameters
% the detected 2D landmarks on the image
mLan = size(lan2d, 1); idxKpts = params.idxKpts;
% camera parameters
phi = params.vecAngle(1); theta = params.vecAngle(2); psi = params.vecAngle(3); 
Rx = [1, 0, 0; 0, cos(phi), -sin(phi); 0, sin(phi), cos(phi)];
Ry = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
Rz = [cos(psi), -sin(psi), 0; sin(psi), cos(psi), 0; 0, 0, 1];
R = Rz*Ry*Rx; 
T = params.T; 
f = params.f; u0 = params.u0; v0 = params.v0;
% shape model parameters
aid = params.mu_id+params.mu_exp; mPts = length(aid)/3; 
Eid = params.pc_id; alpha = params.alpha;
Eexp = params.pc_exp; delta = params.delta;
Mgeo = aid + Eid*alpha + Eexp*delta; Mgeo = reshape(Mgeo, [3, mPts]);
% the histogram of parameters
vec_mParam = params.vec_mParam; mParam = sum(vec_mParam);
mAlpha = vec_mParam(1); mBeta = vec_mParam(2); mDelta = vec_mParam(3); 
mR = vec_mParam(4); mT = vec_mParam(5); mIntrinCam = vec_mParam(6);
m1 = mAlpha+mBeta; m2 = m1+mDelta; m3 = m2+mR; m4 = m3+mT; m5 = m4+mIntrinCam;
% partial derivatives of the rotation matrix with respect to rotation angles
Rx_phi = [0, 0, 0; 0, -sin(phi), -cos(phi); 0, cos(phi), -sin(phi)];
Ry_theta = [-sin(theta), 0, cos(theta); 0, 0, 0; -cos(theta), 0, -sin(theta)];
Rz_psi = [-sin(psi), -cos(psi), 0; cos(psi), -sin(psi), 0; 0, 0, 0];
dR_dphi = Rz*Ry*Rx_phi;
dR_dtheta = Rz*Ry_theta*Rx;
dR_dpsi = Rz_psi*Ry*Rx;

%% Calculate landmark residuals and the Jacobian matrix of residual functions with respect to the feature alignment term
J_rLan = zeros(mLan*2, mParam);
rLan = zeros(mLan*2, 1);
for i = 1:mLan
    idxLan = idxKpts(i);
    Pw = Mgeo(:, idxLan);
    Pc = R*Pw+T;
    Xc = Pc(1); Yc = Pc(2); Zc = Pc(3); Zc2 = Zc^2;  
    idxU = i*2-1; idxV = idxU+1;
    
    %% Feature alignment resiudals
    rLan(idxU) = (f*Xc/Zc+u0)-lan2d(i, 1);
    rLan(idxV) = (f*Yc/Zc+v0)-lan2d(i, 2);
        
    %% Shape parameters
    % partial derivatives with respect to alpha
    idxLanCoor = (idxLan*3-2:idxLan*3);
    dPw_dalpha = Eid(idxLanCoor, 1:mAlpha);
    dPc_dalpha = R*dPw_dalpha;
    J_rLan(idxU, 1:mAlpha) = (dPc_dalpha(1, :)*Zc-Xc*dPc_dalpha(3, :))*f/Zc2; % u
    J_rLan(idxV, 1:mAlpha) = (dPc_dalpha(2, :)*Zc-Yc*dPc_dalpha(3, :))*f/Zc2; % v

    % partial derivatives with respect to delta
    dPw_ddelta = Eexp(idxLanCoor, 1:mDelta);
    dPc_ddelta = R*dPw_ddelta;
    J_rLan(idxU, m1+1:m2) = (dPc_ddelta(1, :)*Zc-Xc*dPc_ddelta(3, :))*f/Zc2; % u
    J_rLan(idxV, m1+1:m2) = (dPc_ddelta(2, :)*Zc-Yc*dPc_ddelta(3, :))*f/Zc2; % v
    
    %% Camera parameters 
    % partial derivatives with respect to rotation angles    
    dPc_dphi = dR_dphi*Pw;
    dPc_dtheta = dR_dtheta*Pw;
    dPc_dpsi = dR_dpsi*Pw;
    dPc_dangle = [dPc_dphi, dPc_dtheta, dPc_dpsi];
    J_rLan(idxU, m2+1:m3) = (dPc_dangle(1, :)*Zc-Xc*dPc_dangle(3, :))*f/Zc2; % u
    J_rLan(idxV, m2+1:m3) = (dPc_dangle(2, :)*Zc-Yc*dPc_dangle(3, :))*f/Zc2; % v

    % partial derivatives with respect to the translation vector     
    J_rLan(idxU, m3+1:m4) = [f/Zc, 0, -(f*Xc)/Zc2]; % u
    J_rLan(idxV, m3+1:m4) = [0, f/Zc, -(f*Yc)/Zc2]; % v
    
    % partial derivatives with respect to the focal length and principal point
    J_rLan(idxU, m4+1:m5) = [Xc/Zc, 1, 0]; % u
    J_rLan(idxV, m4+1:m5) = [Yc/Zc, 0, 1]; % v
end

% assign the normalization factor
normLan_sqrt = sqrt(1/mLan);
J_rLan = normLan_sqrt*J_rLan;
rLan = normLan_sqrt*rLan;

end