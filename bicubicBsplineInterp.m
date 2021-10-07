
function [defIntp] = bicubicBsplineInterp(ImDef,PcoordInt)
% This function is used for bicubic B-spline Iterpolation.
% It is not friendly for the readar beacue I have devoted a lot of effort
% to speed it up. Luckily, it pays back, and is much faster than that using
% for loop.
% Author: Bin Chen;
% E-mail: binchen@kth.se
% Update: 2021-06-04

xInt      = floor(PcoordInt);
deltaX    = PcoordInt-xInt;
numPt     = length(xInt);

MBT       = 1/6*[-1, 3,-3, 1;
                  3,-6, 0, 4;
                 -3, 3, 3, 1;
                  1, 0, 0, 0];
deltaMatX = MBT*[deltaX(1,:).^3;deltaX(1,:).^2;deltaX(1,:);ones(1,numPt)];
deltaMatY = MBT*[deltaX(2,:).^3;deltaX(2,:).^2;deltaX(2,:);ones(1,numPt)];

% the index of calculation pooint in the reference subset
indx      = repmat([xInt(2,:)-2;xInt(2,:)-1;xInt(2,:);xInt(2,:)+1],4,1).*size(ImDef,1)+...
            [repmat(xInt(1,:)-1,4,1);repmat(xInt(1,:),4,1);...
             repmat(xInt(1,:)+1,4,1);repmat(xInt(1,:)+2,4,1)];
D_all      = ImDef(indx);
defIntp   = repmat(deltaMatY,4,1) .* D_all .* [repmat(deltaMatX(1,:),4,1);...
            repmat(deltaMatX(2,:),4,1);repmat(deltaMatX(3,:),4,1);...
            repmat(deltaMatX(4,:),4,1)];
% the interped gray intensity in the deformed image.
defIntp   = sum(defIntp,1)';