function [p,Czncc,Iter,disp]=iterICLM(ImRef,ImDef,pCoord,p,Params)
% IC-LM for first order shape function
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
    localSub = localSub./repmat(ceil(subset./2),length(localSub),1);
    M       = diag([1,1/subset(1),1/subset(2),1,1/subset(1),1/subset(2)]);
else
    M = eye(6);
end

% estimate deltaf of Hessian matrix
nablafx       = gradxImR(pCoord(1)+deltaVecX,pCoord(2)+deltaVecY);
nablafy       = gradyImR(pCoord(1)+deltaVecX,pCoord(2)+deltaVecY);
nablaf        = [nablafx(:),nablafy(:)];
[~,Jacob,H]   = calHessian(localSub,nablaf);
fSubset       = ImRef(pCoord(1)+deltaVecX,pCoord(2)+deltaVecY);
fSubset       = fSubset(:);


deltafVec     = fSubset-mean(fSubset(:));
deltaf        = sqrt(sum(deltafVec.^2));

warP          = Warp(p);
thre          = 1;

% parameters for L-M algorithm
I             = eye(size(H));
lambda        = 1e4;
alpha         = 0.1;
beta          = 10;
Cznssd_0      = 4;


% iteration to estimate p
while thre> 1e-3 && Iter<Params.maxIter || Iter==0
    
    
    gIntep    = warP*localSubHom;
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
        
        Cznssd    = sum((deltafVec/deltaf-deltagVec/deltag).^2);
        
        Iter      = Iter+1;
        if Iter == 1
           lambda = lambda^(Cznssd/Cznssd_0)-1;
        end   
        % update Hessian in L-M
        H_update  = (H+lambda*I);
         % deltap: the increment of p;
        deltap    = -inv(H_update)*Jacob'*delta;
        deltap    = M*deltap;     
        
        
        thre      = sqrt(deltap(1).^2+deltap(4).^2);
        Czncc     = 1-0.5*Cznssd;
        
        % update lambda and warp
        if Cznssd < Cznssd_0
            lambda    = lambda*alpha; % update lambda: damping parameter
            warpdelta = Warp(deltap); % update the warp
            warP      = warP/warpdelta;
            Cznssd_0  = Cznssd;
        else
            lambda = lambda*beta;
        end
        p         = [warP(1,3);warP(1,1)-1;warP(1,2);...
                     warP(2,3);warP(2,1);warP(2,2)-1];
        disp      = p([1,4],:);

    else
        Iter      = Iter+1;
        p         = zeros(6,1);

        disp      = nan(1,2);
        Czncc     = -1;
        break
    end

end

end