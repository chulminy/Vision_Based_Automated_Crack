%% Vision-Based Automated Crack Detection for Bridge Inspection
%
%% Contact
%  Name  : Chul Min Yeum
%  Email : chulminy@gmail.com
%  Please contact me if you have a question or find a bug. You can use a
%  "issues" page in the github.

%% Description
%
% This code is used for the following paper:
% Chul Min Yeum and Shirley J Dyke.“Vision-Based Automated Crack Detection
% for Bridge Inspection.” Computer-Aided Civil and Infrastructure
% Engineering 30, no. 10 (2015): 759–770.
%
% All data and a package code are provided. I did my best to optimize the
% code in terms of speed and readability but, it is still slow to complete
% all steps. For future users for this code, please convert this code to
% C++ or other fater languages.
%
% This code is tested under only window machine so you may find errors in
% other OS, but it may be minor errors related with setting up file pathes
% or toolbox download in "Parameters". I am open to assist you to work this
% code on differt OS.
%
% You can feel free to modify this code but I would recommend cite my
% paper.

%% Reference
% 
% Chul Min Yeum and Shirley J. Dyke. “Vision-based Automated Visual
% Inspection of Large-scale Bridges”. Sixth World Conference on Structural
% Control and Monitoring, July, 2014.
% 
% Chul Min Yeum and Shirley J Dyke.“Vision-Based Automated Crack Detection
% for Bridge Inspection.” Computer-Aided Civil and Infrastructure
% Engineering 30, no. 10 (2015): 759–770.

clear; clc; close all; format shortg; warning off;

%% Setup folders
folderImg           = fullfile(cd(cd('..')),'img');

% training image
folderImgTrOrg      = fullfile(folderImg,'train','org');
folderImgTrResize   = fullfile(folderImg,'train','resize');
folderImgTrData     = fullfile(folderImg,'train','data');

imgFiles  = dir(fullfile(folderImgTrOrg,'*.jpg'));
orgTrImgs = cellfun(@(x) fullfile(folderImgTrOrg,x), ...
                    {imgFiles(:).name},'UniformOutput', false);
rszTrImgs = cellfun(@(x) fullfile(folderImgTrResize,x), ...
                    {imgFiles(:).name},'UniformOutput', false);
clearvars imgFiles;

% test image
folderImgTsOrg      = fullfile(folderImg,'test','org');
folderImgTsResize   = fullfile(folderImg,'test','resize');
folderImgTsData     = fullfile(folderImg,'test','data');
folderImgTsBolt     = fullfile(folderImg,'test','bolt');

imgFiles  = dir(fullfile(folderImgTsOrg,'*.jpg'));
orgTsImgs = cellfun(@(x) fullfile(folderImgTsOrg,x), ...
                    {imgFiles(:).name},'UniformOutput', false);
rszTsImgs = cellfun(@(x) fullfile(folderImgTsResize,x), ...
                    {imgFiles(:).name},'UniformOutput', false);
                
% output image
folderOut   = fullfile(folderImg,'output'); 
folderHome  = fullfile(folderOut,'homepage'); 
clearvars imgFiles folderImg;

%% Include external codes
folder_external = 'external'; mkdir(folder_external);
addpath(genpath(folder_external));

folder_misc = 'misc'; 
addpath(folder_misc);

%% Installation of a VLFEAT toolbox (one-time process)
% If there is error in this code block, you can manually install vlfeat in
% your machine.
% Please check out the following website: http://www.vlfeat.org/
vlfeat_link = 'http://www.vlfeat.org/download/vlfeat-0.9.18-bin.tar.gz';
if ~exist('vl_version','file')
    untar(vlfeat_link,folder_external);
    run(fullfile(folder_external,'vlfeat-0.9.18','toolbox','vl_setup'));
    savepath;
end; clearvars vlfeat_link

%% Installation of piotr's computer vision toolbox
% If there is error in this code block, you can manually install vlfeat in
% your machine. Please check out the following website:
% https://github.com/pdollar/toolbox
pdol_link = 'https://pdollar.github.io/toolbox/archive/piotr_toolbox.zip';
if  ~exist(fullfile(folder_external,'toolbox'),'dir')
    unzip(pdol_link,folder_external);
    addpath(genpath(fullfile(folder_external,'toolbox'))); 
else
    addpath(genpath(fullfile(folder_external,'toolbox'))); 
end; clearvars pdol_link

%% Installation of matlab-util developed by tsogkas
% https://github.com/tsogkas/matlab-utils
util_link = 'https://github.com/tsogkas/matlab-utils/archive/master.zip';
if  ~exist(fullfile(folder_external,'matlab-utils-master'),'dir')
    unzip(util_link,folder_external);
    addpath(genpath(fullfile(folder_external,'matlab-utils-master'))); 
else
    addpath(genpath(fullfile(folder_external,'matlab-utils-master'))); 
end; clearvars util_link

%% Installation of Community Detection Toolbox developed by tsogkas
% https://www.mathworks.com/matlabcentral/fileexchange/45867-community-detection-toolbox/content/ComDetTBv090/Algorithms/GCModulMax3.m
cm_link = 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/45867/versions/1/download/zip';
if  ~exist(fullfile(folder_external,'ComDetTBv090'),'dir')
    unzip(cm_link,folder_external);
    run(fullfile(folder_external,'ComDetTBv090','PathAdd.m')); 
else
    run(fullfile(folder_external,'ComDetTBv090','PathAdd.m')); 
end; clearvars folder_external cm_link

%% Parameters

% # of images
nTsImg = numel(orgTsImgs);

% img resizing scale (resize)
imgRszScale = 0.4;

% training window patch size
nHaar           = 500; % # of Haar patches
patchSize       = 64;  % a path size for converted images
patchSizeOrg    = patchSize/imgRszScale;
patchStride     = 8;   % stride windows

% # of channel features 
nChannel     = 11;
% nLUVChn      = 2;           % LUV and remove L component
% nHSVChn      = 2;           % HSV and remove V component
% nOrients     = 6;           % histogram orientation. 
% nGradHistChn = nOrients;    
% nGradImgChn  = 1;           % gradient image
% nChannel     = nLUVChn + nHSVChn + nGradHistChn + nGradImgChn;

% # of features
nFeature    = nChannel * nHaar;

% Adaboost learning round (# of wek classifier)
nRounds = 200;

% gaussian filter (remove high frequency noise initially)
gaussF = fspecial('gaussian', 3, 0.5);

% image downsacle
scaleFactor = 1/2^(1/5);
nScale      = 6;
imgScale    = scaleFactor.^(0:nScale-1);

% negative image generation
nNegMul         = 5;   % nNeg = nNegMul * nPos 

% select negative window size among 'negWinSize'
negWinSize   = patchSizeOrg:4/imgRszScale:patchSizeOrg/imgScale(end);

% scale multiplication (window size multiplication)
scaleMul     = 1.5;

% IOU threshold (for merging detections)
thresh_IOU   = 0.4;

% # of match images (for matching objects: next numCompImg images are
% corresponded based on the assumption that the images are captured
% consecutively.
numCompImg   = 5;

% # of match images (for grouping)
nMinFeatMatch = 3; 
nMinPatch     = 5;

% Parameters for crack detection ------------------------------------------
% Frangi parameters
Param.FrangiScaleRange    = [1 3];
Param.FrangiScaleRatio    = 0.2;
Param.FrangiBetaOne       = 0.5;
Param.FrangiBetaTwo       = 0;
Param.verbose             = false;
Param.BlackWhite          = false;

% distance mask 
Param.borderMargin        = 6;
Param.minRatio            = 0.5;

% crack supurious edge
Param.minEdgeArea         = 30;
Param.FrangiFilterThresh  = 0.85;

% theta for radon transformation
Param.theta               = 0:180;

% edge shape descriptor
Param.rgd_ratioMajMin     = 8;        % regional descriptor
Param.rgd_minorAxis       = 10;