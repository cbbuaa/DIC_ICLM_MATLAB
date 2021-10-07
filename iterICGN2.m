function [p,Czncc,Iter,disp] = iterICGN2(ImRef,ImDef,pCoord,p,Params)
% IC-GN for 2nd order shape function
% Author: Bin Chen;
% E-mail: binchen@kth.se
% Update: 2021-06-04

Iter          = 0;
subset        = Params.subset;
localSubHom   = Params.localSubHom;
deltaVecX     = Params.deltaVecX;
deltaVecY     = Params.deltaVecY;
gradxImR      = Params.gradxImR;
gradyImR      = Params.gradyImR;
localSub      = Params.localSub;
numPt         = prod(subset);


% normalization
if Params.Normalization
    localSub  = localSub./repmat(ceil(subset./2),length(localSub),1);
    M         = diag([1,1/subset(1),1/subset(2),1/subset(1)^2,1/prod(subset),1/subset(2)^2,...
             1,1/subset(1),1/subset(2),1/subset(1)^2,1/prod(subset),1/subset(2)^2]);
else
    M         = eye(12);
end

% estimate deltaf of Hessian matrix
nablafx       = gradxImR(pCoord(1)+deltaVecX,pCoord(2)+deltaVecY);
nablafy       = gradyImR(pCoord(1)+deltaVecX,pCoord(2)+deltaVecY);

nablaf        = [nablafx(:),nablafy(:)];
[invH,Jacob]  = calHessian2(localSub,nablaf,numPt);

%% deltaf
fSubset       = ImRef(pCoord(1)+deltaVecX,pCoord(2)+deltaVecY);
fSubset       = fSubset(:);

deltafVec     = fSubset-mean(fSubset(:));
deltaf        = sqrt(sum(deltafVec.^2));
invHJacob     = invH*Jacob';

warP          = Warp(p);
thre          = 1;

% iteration to estimate p
while thre>1e-3 && Iter< Params.maxIter || Iter==0
    % OK
    gIntep    = warP([4,5,6],:)*[localSubHom(1,:).^2;localSubHom(1,:).*localSubHom(2,:);localSubHom(2,:).^2;localSubHom];
    PcoordInt = pCoord + gIntep-[0;0;1];
    if prod(min(PcoordInt(1:2,:),[],2)>[3;3]) && ...
            prod(max(PcoordInt(1:2,:),[],2)<[Params.sizeX-3;Params.sizeY-3])
%       B-spline interpolation, the following are two choice:
%       1. based on Matlab, and the code has been vectorized for
%       effeciency;
%       2. based on C++ mex file, is aroud 4 times faster.
%         defIntp   = bicubicBsplineInterp(ImDef,PcoordInt)
        defIntp   = BicubicBsplineInterp(ImDef,PcoordInt);        
        
        deltagVec = defIntp-mean(defIntp);
        deltag    = sqrt(sum(deltagVec.^2));
        delta     = deltafVec-deltaf/deltag*deltagVec;
        deltap    = -invHJacob*delta;
        
        deltap    = M*deltap;

        Iter      = Iter+1;
  
        thre      = sqrt(deltap(1).^2+deltap(7).^2);
        
        Cznssd    = sum((deltafVec/deltaf-deltagVec/deltag).^2);
        Czncc     = 1-0.5*Cznssd;
        
        warpdelta = Warp(deltap);
        warP      = warP/warpdelta;
        
        p         = [warP(4,6);warP(4,4)-1;warP(4,5);2*warP(4,1);warP(4,2);2*warP(4,3);...
                     warP(5,6);warP(5,4);warP(5,5)-1;2*warP(5,1);warP(5,2);2*warP(5,3)];

        disp      = p([1,7],:);

    else
        Iter      = Iter+1;
        p         = zeros(12,1);

        disp      = nan(1,2);
        Cznssd    = sum((deltafVec/deltaf-deltagVec/deltag).^2)/numPt;
        Czncc     = -1;
        break
    end
end