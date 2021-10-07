function [Params]=calcuPt(Params)
% this function is used for generate a set of calculation points
% and show ROI and subset on the figure
% Author: Bin Chen;
% E-mail: binchen@kth.se
% Update: 2021-06-05

% parameters
half_subset  = Params.half_subset;
Step         = Params.Step;

%% get two corner points of ROI, and plot a rectangle
hold on,

if Params.fixed_seedPts
    x = [100,Params.sizeX-100];
    y = [100,Params.sizeY-100];
else
    [y(1),x(1)]  = ginput(1);
    plot(y(1),x(1),'+r','MarkerSize',8);
    hold on
    [y(2),x(2)]  = ginput(1);
    hold on
    plot(y(2),x(2),'+r','MarkerSize',8);
end
x            = round(x);
y            = round(y);
xROI         = sort(x);
yROI         = sort(y);

plot([yROI(1),yROI(1)],[xROI(1),xROI(2)],'-r','LineWidth',1);
plot([yROI(2),yROI(2)],[xROI(1),xROI(2)],'-r','LineWidth',1);
plot([yROI(1),yROI(2)],[xROI(1),xROI(1)],'-r','LineWidth',1);
plot([yROI(1),yROI(2)],[xROI(2),xROI(2)],'-r','LineWidth',1);

%% generate a set of calculaiton points with space equal to "step"
calPtX          = (xROI(1)+half_subset(1)):Step:(xROI(2)-half_subset(1));
calPtY          = (yROI(1)+half_subset(2)):Step:(yROI(2)-half_subset(2));
[PMeshX,PMeshY] = ndgrid(calPtX,calPtY);

% number of calculation points
Lx              = length(calPtX); % point number along x direction
Ly              = length(calPtY); % point number along y direction
[refX,refY]     = ndgrid(calPtX,calPtY); % grid of calulation point
comptPoints     = [refX(:),refY(:)]; % the calculation points

%% the initial point
if Params.fixed_seedPts
    % the central point is select as the initial point
    x               = Params.sizeX./2;
    y               = Params.sizeY./2;
    dist            = sqrt((PMeshX(:)-x).^2+(PMeshY(:)-y).^2);
    [~,indx]        = min(dist);
    InitP           = [PMeshX(indx(1));PMeshY(indx(1));1];
else
    % manually select an initial point
    [x,y]           = ginput(1);
    dist            = sqrt((PMeshX(:)-y).^2+(PMeshY(:)-x).^2);
    [~,indx]        = min(dist);
    InitP           = [PMeshX(indx(1));PMeshY(indx(1));1];
end
% plot it on the image
plot(InitP(2),InitP(1),'r*','MarkerSize',12);
drawnow;
hold on

% save all points into structure: Params.
Params.Lx          = Lx;
Params.Ly          = Ly;
Params.InitP       = InitP;
Params.IndxInitP   = indx;
Params.calPtX      = calPtX;
Params.calPtY      = calPtY;
Params.NumCalPt    = Lx * Ly;
Params.comptPoints = comptPoints;
