function [imgCrop, lan2dCrop] = cropImg(img, lan2d)
% Crop the face image according to the landmarks

wImg = size(img, 2);
hImg = size(img, 1);

minXY = min(lan2d);
minX = minXY(1);minY = minXY(2);
maxXY = max(lan2d);
maxX = maxXY(1);maxY = maxXY(2);
wLan2d = maxX - minX;
hLan2d = maxY - minY;

minXCrop = round(max(minX - wLan2d/3, 0));
minYCrop = round(max(minY - hLan2d/3, 0));
wImgCrop = round(min(wLan2d/3*5, wImg-minXCrop-1));
hImgCrop = round(min(hLan2d/3*5, hImg-minYCrop-1));

if wImgCrop < 5 || hImgCrop < 5 || wImgCrop*hImgCrop < 100
    imgCrop = []; % the face is too small
    lan2dCrop = [];
else
    % crop 
    imgCrop = imcrop(img, [minXCrop, minYCrop, wImgCrop, hImgCrop]);

    % scale
    s = max(1, sqrt(wImgCrop*hImgCrop/170000));
    if s > 1
        imgCrop = imresize(imgCrop, 1/s, 'nearest');
    end
    nLan = size(lan2d, 1);
    lan2dCrop(:, 1) = (lan2d(:, 1)-minXCrop*ones(nLan, 1))/s;
    lan2dCrop(:, 2) = (lan2d(:, 2)-minYCrop*ones(nLan, 1))/s;
end

end



	