function strain = strainEstPLS(Disp,Params)
% Calculate the strain fields based on the displacement based on pointwise
% least-squares.

% Ref: B. Pan, Full-field strain measurement using a two-dimensional
% Savitzky-Golay digital differentiator in digital image correlation,
% Opt. Eng. 46 (2007) 033601.

% Author: Bin Chen;
% E-mail: binchen@kth.se
% Update: 2021-06-04

strainWinSize = Params.strainWin;
X             = Params.comptPoints;

ROIsize       = [Params.Lx,Params.Ly];
A             = NaN(prod(ROIsize),3);
B             = A;
halfWsize     = floor(strainWinSize/2);

x = reshape(X(:,1),ROIsize);    y = reshape(X(:,2),ROIsize);
u = reshape(Disp(:,1),ROIsize); v = reshape(Disp(:,2),ROIsize);
for i = 1 : ROIsize(1)
    for j = 1 : ROIsize(2)

        xCen                = x(i,j);
        yCen                = y(i,j);


        xVec            = max(1,i-halfWsize) : min(ROIsize(1),i+halfWsize);
        yVec            = max(1,j-halfWsize) : min(ROIsize(2),j+halfWsize);

        uWin            = u(xVec,yVec);
        vWin            = v(xVec,yVec);
        %             wWin            = w(xVec,yVec);
        xWin            = x(xVec,yVec);
        yWin            = y(xVec,yVec);

        deltax          = xWin(:)-xCen;
        deltay          = yWin(:)-yCen;

        indNan          = find(isnan(deltax));
        uWin(indNan)    = [];
        vWin(indNan)    = [];
        deltax(indNan)  = [];
        deltay(indNan)  = [];
        
        
%         COSNT           = 2000;
%         coeff           = [ones(numel(uWin),1),deltax+COSNT,deltay+COSNT]';        
%         A((j-1)*ROIsize(1)+i,:) = robustfit(coeff(2:3,:)',uWin','welsch');
%         B((j-1)*ROIsize(1)+i,:) = robustfit(coeff(2:3,:)',vWin','welsch');

%         COSNT           = 2000;
%         coeff           = [ones(numel(uWin),1),deltax+COSNT,deltay+COSNT]';
%         coeffMat        = (coeff*coeff')^-1;
%         A((j-1)*ROIsize(1)+i,:) = coeffMat*(coeff*uWin');
%         B((j-1)*ROIsize(1)+i,:) = coeffMat*(coeff*vWin');
%         
        coeff           = [ones(numel(uWin),1),deltax,deltay]';
        coeffMat        = (coeff*coeff')^-1;
        A((j-1)*ROIsize(1)+i,:) = coeffMat*(coeff*uWin');
        B((j-1)*ROIsize(1)+i,:) = coeffMat*(coeff*vWin');
    end
end




% output the strain
exx      = A(:,2)*1e6;
eyy      = B(:,3)*1e6;
exy      = (A(:,3)+B(:,2))*1e6;
strain   = [exx,eyy,exy];