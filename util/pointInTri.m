function [isInside, a, b, c] = pointInTri(pt2d, matPt2dTri)
% Test if a point locate inside a triangle based on barycentric coodinate
% system

x = pt2d(1); y = pt2d(2);
x1 = matPt2dTri(1, 1); y1 = matPt2dTri(1, 2);
x2 = matPt2dTri(2, 1); y2 = matPt2dTri(2, 2);
x3 = matPt2dTri(3, 1); y3 = matPt2dTri(3, 2);
denominator = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3);
a = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / denominator;
b = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / denominator;
c = 1 - a - b;

if 0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c<= 1
    isInside = true;
else
    isInside = false;
end

end

