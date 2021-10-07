% Before run it, please do remember to set the parameters in 'paramset.m'
% NOTE: x is along vertical, and y is along horizontal
% Author: Bin Chen;
% E-mail: binchen@kth.se
% Update: 2021-06-04

clear; close all; clc;

% first select the refenrence image, and deformed image
[refFile, imgPath]   = uigetfile({'*.bmp;*.tiff;*.tif'},'Select reference image');
[dataFileAllTemp, ~] = uigetfile({'*.bmp;*.tiff;*.tif'},'MultiSelect','on',...
    'Select deformed images');
if ~iscell(dataFileAllTemp)
    dataFileAll{1}   = dataFileAllTemp;
else
    dataFileAll      = dataFileAllTemp;
end

% Set parameters
Params               = paramset();
   
% Read and show reference image
subplot(1,2,1);
title(refFile);
ImRef0               = double(imread(fullfile(imgPath,refFile)));
h                    = fspecial('gaussian',5,1);
ImRef                = imfilter(ImRef0,h);
imshow(repmat(uint8(ImRef),1,1,3));
% normalize the image
if Params.Normalization
ImRef                = ImRef - mean(ImRef(:));
ImRef                = ImRef./max(abs(ImRef(:)));
end
Params.sizeX         = size(ImRef,1);
Params.sizeY         = size(ImRef,2);

% Select calculation points on the reference image
Params               = calcuPt(Params);
 % number of calculation points

% Gradient of reference image, needs noly be calculated once=
[gradxImR, gradyImR] = gradImg(ImRef);
Params.gradxImR      = gradxImR;
Params.gradyImR      = gradyImR;

% if the data files have been generated, you can skip the matching
% by setting reRun = 0
reRun = 1;
% by changing "for" to "parfor" in this loop, it can make full use of all
% cores of your computer.
for i = 1:length(dataFileAll)
    defFile            = dataFileAll{i};
    [~,dataFileName,~] = fileparts(defFile);
    dataFile           = fullfile(imgPath,[dataFileName,'.mat']);
    if reRun || ~exist(dataFile,'file')
        fprintf('Target image: %s;\n',defFile);
        % Matching
        [disp,strain,ZNCC,iterNum,outParams] = DIC_main(imgPath,defFile,ImRef,Params);
        parsave(dataFile,disp,strain,ZNCC,iterNum,outParams); 
        % select the last input among {u, v, exx, eyy}
        plotOnImg(outParams, ImRef0, outParams.calPtX, outParams.calPtY, disp, strain, 'exx')  
    else
        fprintf('Image %s has been registered;\n',defFile); 
    end                           
end

