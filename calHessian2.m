function [invH,Jacobian,H] = calHessian2(localSub,nablaf,numPt)
% Estimate the Hessian matrix for first order shape function;
% Author: Bin Chen;
% E-mail: binchen@kth.se
% Update: 2021-06-05

deltaW2P    = [ones(numPt,1),localSub(:,1),localSub(:,2),1/2*localSub(:,1).^2,...
               localSub(:,1).*localSub(:,2),1/2*localSub(:,2).^2];
Jacobian    = [repmat(nablaf(:,1),1,6).*deltaW2P,repmat(nablaf(:,2),1,6).*deltaW2P];
H           = Jacobian'*Jacobian;
invH        = inv(H);