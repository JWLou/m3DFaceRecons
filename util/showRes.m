function showRes(params, img, lan2d)
% Show 3D shape and textures with current parameters
clf;

R = rotationMatrix(params.vecAngle(1), params.vecAngle(2), params.vecAngle(3));
T = params.T;
face3dMean = params.mu_id + params.mu_exp;
face3dEst = face3dMean + params.pc_id*params.alpha + params.pc_exp*params.delta;
nPts = length(face3dMean)/3;
face3dEst = reshape(face3dEst, [3, nPts]);
face3dEstR = R*face3dEst+repmat(T, [1, nPts]);
lan3dEstR = face3dEstR(:, params.idxKpts);

subplot(2,2,1);
imshow(img);
title('Input face image');
hold on;
% plot(lan2d(:, 1), lan2d(:, 2), 'g.');
hold off;
subplot(2,2,2);
fprintf('Synthesizing the face image...\n');
[imgSyn, mapTriPix] = synFaceImg(params, img, lan2d);
imgSyn = uint8(imgSyn);
imshow(imgSyn);
title('Synthesized face image')
subplot(2,2,3);
face3dMean = reshape(face3dMean, [3, nPts]);
face3dMeanR = R*face3dMean+repmat(T, [1, nPts]);
land3dMeanR = face3dMeanR(:, params.idxKpts);
[~, ~, idxTriPix] = find(mapTriPix);
drawHead(face3dMeanR, params.tri(:, params.idxFaceTri), land3dMeanR);
title('Initial mean 3D face');
% trimesh(params.tri', face3dMean(1, :), face3dMean(2, :), face3dMean(3, :));
subplot(2,2,4);
% trimesh(params.tri(:, idxTriPix)', face3dEst(1, :), face3dEst(2, :), face3dEst(3, :));
drawHead(face3dEstR, params.tri(:, params.idxFaceTri), lan3dEstR);
title('Recovered 3D face');
pause;

end
