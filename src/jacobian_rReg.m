function [J_rReg, rReg] = jacobian_rReg(params)
% Calculate the regularization residuals with respect to the current parameters and the Jacobian matrix 
% of regularization residuals with respect to the parameters based on IRLS solver - [alpha, beta, delta, 
% vecAngle(phi(angleX), theta(angleY), psi(angleZ)), T, f, (u0, v0), gamma(3*9)]

%% Parameters
% shape and reflectance parameters' standard deviations
sigmaId = params.ev_id; sigmaExp = params.ev_exp; sigmaAlb = params.ev_alb; 
alpha = params.alpha; delta = params.delta; beta = params.beta;
% the histogram of parameters
vec_mParam = params.vec_mParam; mParam = sum(vec_mParam);
mAlpha = vec_mParam(1); mBeta = vec_mParam(2); mDelta = vec_mParam(3);
mReg = mAlpha+mBeta+mDelta;

%% Calculate the regularization residuals with respect to the current parameters and the Jacobian matrix of residual functions with respect to the statistical regularization term
J_rReg = zeros(mReg, mParam);
rReg = zeros(mReg, 1);
for i = 1:mReg
    if i <= mAlpha
        idxTmp = i;
        rReg(i) = alpha(idxTmp)/sigmaId(idxTmp);
        
       %% the partial derivative with respect to identity weights
        J_rReg(i, i) = 1/sigmaId(idxTmp);    
    elseif i <= mAlpha+mBeta
        idxTmp = i-mAlpha;
        rReg(i) = beta(idxTmp)/sigmaAlb(idxTmp);
        
       %% the partial derivative with respect to albedo weights
        J_rReg(i, i) = 1/sigmaAlb(idxTmp);
    else
        idxTmp = i-mAlpha-mBeta;
        rReg(i) = delta(idxTmp)/sigmaExp(idxTmp);
        
       %% the partial derivative with respect to expression weights
        J_rReg(i, i) = 1/sigmaExp(idxTmp);
    end
end

end




