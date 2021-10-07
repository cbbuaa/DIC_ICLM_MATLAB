function Params = paramset(Params)
% This function is used for setting all parameters
% Author: Bin Chen;
% E-mail: binchen@kth.se
% Update: 2021-06-04

%% -------------Check these parameters before usage-------------%%
if nargin == 0
    % set the subset size;
    subset           = 31*[1,1];
else
    subset           = Params.subset;
end
Step                 = 5;
strainWin            = [9,9];

% IterMethod: choose different iteration methods, including 'IC-GN','IC-LM',
% 'IC-GN2','IC-LM2','IC-GN-RG','IC-LM-RG','IC-GN2-RG','IC-LM2-RG'.
IterMethod           = 'IC-LM2-RG';
Normalization        = 1; % 1 with normalization, 0 without normalization
maxIter              = 15; % maximum iteration number
thre                 = 1e-3; % theroshold of iteration
fixed_seedPts        = 0; % 0: manually select; 1: fixed to the center point
% InitMethod: 0 denotes no initial value is set,1 denotes corse-to-fine, 
% 2 denotes Fast Fourier Transformationo, 3 denotes manually selection
InitMethod           = 2;

%% -------------------------------------------------------------- %%
half_subset          = floor(subset/2);
% local coordinate of a subset
deltaVecX            = -half_subset(1) : half_subset(1);
deltaVecY            = -half_subset(2) : half_subset(2);
[deltax,deltay]      = ndgrid(deltaVecX,deltaVecY);
localSubHom          = [deltax(:)';deltay(:)';ones(1,prod(subset))];

% Local coordinate in subset
localSub             = [localSubHom(1,:)',localSubHom(2,:)'];
% Print the key parameters on the screen.
fprintf('Matching method: %s;\n', IterMethod);
fprintf('Subset size: %d * %d;\n', subset(1),subset(2));
fprintf('Strain window size: %d * %d;\n', strainWin(1),strainWin(2));
fprintf('Step: %d;\n\n', Step);


%% Save all results into a structure
Params.fixed_seedPts = fixed_seedPts;
Params.thre          = thre;
Params.subset        = subset;
Params.Step          = Step;
Params.maxIter       = maxIter;
Params.half_subset   = half_subset;
Params.localSubHom   = localSubHom;
Params.deltaVecX     = deltaVecX;
Params.deltaVecY     = deltaVecY;
Params.strainWin     = strainWin;
Params.InitMethod    = InitMethod;
Params.IterMethod    = IterMethod;
Params.Normalization = Normalization;
Params.localSub      = localSub;

