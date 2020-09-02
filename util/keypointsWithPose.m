function idxKptsPose = keypointsWithPose(angleZ, angleY, angleX, shape3d, cheekContours, idxKpts)
% Update cheek keypoint indices across poses

R = angle2dcm(angleZ, angleY, angleX);
projectedShape3d = R * shape3d;
idxKptsPose = idxKpts;
    
% Get the outermost point on the contour line
if(angleY >= 0)
    for i = 1:8
        [~, min_idx] = min(projectedShape3d(1, cheekContours{i}));
        idxKptsPose(i) = cheekContours{i}(min_idx);
    end
else
    for i = 10:17
        [~, max_idx] = max(projectedShape3d(1, cheekContours{i}));
        idxKptsPose(i) = cheekContours{i}(max_idx);
    end
end
    
end