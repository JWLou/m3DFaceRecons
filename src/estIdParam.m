function params = estIdParam(params, lan2d)
% Estimate identity and expression coefficients with a known camera model
% and the perspective-n-points constraint using Gauss-Newton algorithm

%% Parameters
% the detected 2D landmarks on the image
mLan = size(lan2d, 1); idxKpts = params.idxKpts;
% camera parameters
phi = params.vecAngle(1); theta = params.vecAngle(2); psi = params.vecAngle(3); 
Rx = [1, 0, 0; 0, cos(phi), -sin(phi); 0, sin(phi), cos(phi)];
Ry = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
Rz = [cos(psi), -sin(psi), 0; sin(psi), cos(psi), 0; 0, 0, 1];
R = Rz*Ry*Rx; 
T = params.T; f = params.f; u0 = params.u0; v0 = params.v0;
% shape model parameters
aid = params.mu_id+params.mu_exp; mPts = length(aid)/3; 
Eid = params.pc_id; alpha = params.alpha; sigmaId = params.ev_id; 
Eexp = params.pc_exp; delta = params.delta;
Mgeo = aid + Eid*alpha + Eexp*delta; Mgeo = reshape(Mgeo, [3, mPts]);
% the histogram of parameters
vec_mParam = params.vec_mParam;
mAlpha = vec_mParam(1); 

%% Gauss-Newton step
wLan = 100; wLan_sqrt = sqrt(wLan);
wReg = 2.5*10^1; wReg_sqrt = sqrt(wReg);
J_rLan = zeros(mLan*2, mAlpha);
rLan = zeros(mLan*2, 1);
J_rReg = zeros(mAlpha, mAlpha);
rReg = zeros(mAlpha, 1);

%% Landmark term
for i = 1:mLan
    idxLan = idxKpts(i);
    Pw = Mgeo(:, idxLan);
    Pc = R*Pw+T;
    Xc = Pc(1); Yc = Pc(2); Zc = Pc(3); Zc2 = Zc^2;  
    idxU = i*2-1; idxV = idxU+1;

    % feature alignment resiudals
    rLan(idxU) = (f*Xc/Zc+u0)-lan2d(i, 1);
    rLan(idxV) = (f*Yc/Zc+v0)-lan2d(i, 2);

    % Shape parameters
    % partial derivatives with respect to alpha
    idxLanCoor = (idxLan*3-2:idxLan*3);
    dPw_dalpha = Eid(idxLanCoor, 1:mAlpha);
    dPc_dalpha = R*dPw_dalpha;
    J_rLan(idxU, 1:mAlpha) = (dPc_dalpha(1, :)*Zc-Xc*dPc_dalpha(3, :))*f/Zc2; % u
    J_rLan(idxV, 1:mAlpha) = (dPc_dalpha(2, :)*Zc-Yc*dPc_dalpha(3, :))*f/Zc2; % v

end
% assign the normalization factor
normLan_sqrt = sqrt(1/mLan);
J_rLan = normLan_sqrt*J_rLan;
rLan = normLan_sqrt*rLan;

% Regularization term
for i = 1:mAlpha
    rReg(i) = alpha(i)/sigmaId(i);
    % the partial derivative with respect to identity weights
    J_rReg(i, i) = 1/sigmaId(i);        
end

 % assign the term weight
 J_rLan = wLan_sqrt*J_rLan; rLan = wLan_sqrt*rLan;
 J_rReg = wReg_sqrt*J_rReg; rReg = wReg_sqrt*rReg;

% the overall Jacobian and residuals
 J_rAll = [J_rLan; J_rReg];
 rAll = [rLan; rReg];

 % calculate the analytical solution
 A = J_rAll'*J_rAll; b = -J_rAll'*rAll;
 deltaParam = zeros(mAlpha, 1);
 deltaParam = solverPCG(A, b, deltaParam, params.nIterPCG);

 params.alpha = alpha+deltaParam;


end