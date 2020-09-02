function params = initParam(img)
% Initialize all the parameters involved in the optimization process

fprintf('Initializing the parametric face model...\n');

modelId = load('./model/modelId.mat');
modelExp = load('./model/modelExp.mat');
modelTex = load('./model/modelTex.mat');
load('./model/cheekContours.mat');

%% Energy term weights
params.wColor = 1;
params.wLan = 10; 
params.wReg = 2.5*10^-5;

%% Shape model parameters
nPCid = 80; % identity
params.pc_id = double(modelId.pc_id(:, 1:nPCid));
params.mu_id = double(modelId.mu_id);
params.ev_id = double(modelId.ev_id(1:nPCid));
params.alpha =  zeros(nPCid, 1);
nPCexp = 29; % expression
params.mu_exp = double(modelExp.mu_exp); 
params.pc_exp = double(modelExp.pc_exp(:, 1:nPCexp));
params.ev_exp = double(modelExp.ev_exp(1:nPCexp));
params.delta = zeros(nPCexp, 1);
params.tri = modelId.tri; % triangle
params.idxFaceTri = modelId.idxFaceTri2; 

%% Reflectance(albedo) model parameters
nPCtex = 80;
params.mu_tex = double(modelTex.mu_tex); 
params.pc_tex = double(modelTex.pc_tex(:, 1:nPCtex));
params.ev_alb = double(modelTex.ev_tex(1:nPCtex)); 
params.beta = zeros(nPCtex, 1); 

%% Illumination model parameters
params.gamma = ones(27, 1); 

%% Camera parameters
imgWidth = size(img, 2); 
imgHeight = size(img, 1);
params.f = 1000; % intrinsic
params.u0 = imgWidth/2;
params.v0 = imgHeight/2;
params.vecAngle = [0; 0; 0]; % extrinsic ;
params.T = [0; 0; 0];

%% The histogram of parameters
mAlpha = nPCid; mBeta = nPCtex; mDelta = nPCexp; 
mR = 3; mT = 3; mIntrinCam = 3; mGamma = 27;
params.vec_mParam = [mAlpha, mBeta, mDelta, mR, mT, mIntrinCam, mGamma]; 

%% IRLS and PCG solver iterations
params.nIterIRLS = 7;
params.nIterPCG = 4;

%% Keypoints 
params.idxKpts = modelId.idxKpts;
params.cheekContours = cheekContours;

clear modelId modelExp modelTex cheekContours;

end