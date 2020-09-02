function params = main3DFaceRecons
% Reconstruct 3D face from a single image

clear;
fprintf('Reconstruction begins...\n');

%% Initialization
img = imread('.\data\test.png');
load('.\data\test_pts.mat'); lan2d = pts;
[img, lan2d] = cropImg(img, lan2d);

% initialize parameters
params = initParam(img); 
params = estCamParam(params, lan2d);
% params = estGeoParam(params, lan2d);
% params = estExpParam(params, lan2d);
params = estIllAlbParam(params, img);

% Gaussian pyramid with 3 levels
imgR1 = impyramid(img, 'reduce'); lan2dR1 = 0.5*lan2d;
imgR2 = impyramid(imgR1, 'reduce'); lan2dR2 = 0.5*lan2dR1;

%% Call the IRLS solve
% image pyramid - level 2
fprintf('\nImage pyramid - level 2\n');
params.f = params.f/4; 
params.u0 = params.u0/4; params.v0 = params.v0/4; 
params.nIterIRLS = 7; 
params = solverIRLS(params, imgR2, lan2dR2);

% image pyramid - level 1
fprintf('\nImage pyramid - level 1\n');
params.f = 2*params.f; 
params.u0 = 2*params.u0; params.v0 = 2*params.v0;
params.nIterIRLS = 5;
params = solverIRLS(params, imgR1, lan2dR1);

% image pyramid - level 0
fprintf('\nImage pyramid - level 0\n');
params.f = 2*params.f; 
params.u0 = 2*params.u0; params.v0 = 2*params.v0;
params.nIterIRLS = 3;
params = solverIRLS(params, img, lan2d);

save3d(params, 'test.ply');

end

function save3d(params, fpath)
% Save the recontructed and textured face

verts = params.mu_id + params.mu_exp + params.pc_id*params.alpha + params.pc_exp*params.delta;
verts = reshape(verts, 3, length(verts)/3);
% tex = params.mu_tex + params.pc_tex*params.beta;
% tex = reshape(tex, 3, length(tex)/3);
tex = calVertTex(params);
tex(tex > 255) = 255;
tex(tex < 0) = 0;
tri = params.tri;

plywrite(fpath, tri', verts', uint8(tex'));

end