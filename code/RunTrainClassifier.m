%% Description 
%
% This file is to train a classifier for bolt patch dection. First, using
% on positive and negative patches annotated from "RunLabeling.m", a bolt
% patch classifier is trained. Next, you will apply this classifier to the
% entire areas of all training images and detect positive and negative
% patches. Hard negative patches means that the patches classified as
% negative near the objects but they are true bolt patches. All hard
% negative patches are included in "pos" samples and similar number of
% negative patches are included in "neg" samples for training a new and
% better bolt classifier. Intersect-over-Union (IOU) is used for detect
% these positive and negative samples. 
%
% To increase computational efficiency, I tweaked some codes, which may
% not be intuitively understood. If you don't understand the code, please
% email me and I will make a better description next to the corresponding
% code lines.
%
Parameters;

%% Generate training patches
% The original training images are resized for rapid training.
if ~exist(fullfile(folderImgTrData,'DataTrGTPatch.mat'),'file')
    
    % Users must annotate bolts on the original images prior to run this
    % code(see RunLabeling.m)
    load(fullfile(folderImgTrData,'DataLabelAll.mat'),'boltTrLabel');
    nTrImg = numel(boltTrLabel); % # of training images
    
    count = 1;
    for ii=1:nTrImg
        if ~exist(rszTrImgs{ii},'file')
            img = imresize(imread(orgTrImgs{ii}),imgRszScale);
            img = imfilter(img, gaussF); % remove noise
            imwrite(img,rszTrImgs{ii});
        else
            img = imread(rszTrImgs{ii});
        end
        clssLabel   = boltTrLabel(ii).clssLabel;
        rszBbox     = (boltTrLabel(ii).bbox)*imgRszScale;
        
        nbbox       = size(rszBbox,1);
        for jj=1:nbbox
            points = bbox2points(rszBbox(jj,:)); 
            % please refer the following:
            % https://www.mathworks.com/help/vision/ref/bbox2points.html
            
            % remove the patch overlapping with image boundary
            if  any(points(1,:) < 1) || ...
                    (points(4,1) > size(img,2)) || ...
                    (points(4,2) > size(img,1))
                continue;
            end
            
            boltGTPatch(count).trImgIdx       = ii; % training image index
            boltGTPatch(count).rszOrgBbox     = rszBbox(jj,:);
            boltGTPatch(count).patch          = imresize(imcrop(img,...
                rszBbox(jj,:)), [patchSize patchSize]); % bolt GT patch
            boltGTPatch(count).clss           = clssLabel(jj);
            
            count = count + 1;
        end
    end
    save(fullfile(folderImgTrData,'DataTrGTPatch.mat'),'boltGTPatch');
else
    load(fullfile(folderImgTrData,'DataTrGTPatch.mat'),'boltGTPatch');
end

%% Generate Haar patches
if ~exist(fullfile(folderImgTrData,'DataHaar.mat'),'file');
    haarPatch = cell(nHaar,1);
    for ii= 1:nHaar
        [posMat, operMat]   = CompHaarPatch([1 1], patchSize, patchSize);
        haarPatch{ii,1}     = [posMat operMat];
        
        % Indexes allocation for fast computation haar-like feature using
        % on a summed area table
        % (https://en.wikipedia.org/wiki/Haar-like_features)
        haarComp  = zeros(patchSize+1,patchSize+1);
        haarComp(sub2ind(size(haarComp),posMat(:,2),posMat(:,1)))= operMat;
        
        haarPatch{ii,2} = haarComp;
    end; clearvars haarComp ii jj posMat operMat;
    tmp = repmat(1:nChannel,nHaar,1); featIdx.channel  = tmp(:)';
    featIdx.haar     = repmat(1:nHaar,1,nChannel);
    
    % for example: a feature is computed using a featIdx.haar(i) haar
    % patch on featIdx.channel(i) channel image
    
    save(fullfile(folderImgTrData,'DataHaar.mat'),'haarPatch','featIdx');
else
    load(fullfile(folderImgTrData,'DataHaar.mat'),'haarPatch','featIdx');
end

%% Compute features using labeled patches
if ~exist(fullfile(folderImgTrData,'DataTrPatchFeat.mat'),'file');
    nTrPatch        = numel(boltGTPatch);
    trFeatMat       = zeros(nTrPatch,nFeature);
    trChnPatchIdx   = zeros(nTrPatch,1);
    trHaarPatchIdx  = zeros(nTrPatch,1);
    
    for ii=1:nTrPatch
        % compute integral channel images
        itgPatchSet      = CompIntChnImg(boltGTPatch(ii).patch);
        
        % compute features from integral channel images
        trFeatMat(ii,:)  = CompFeatChnPatch(itgPatchSet,haarPatch, nHaar);
    end
    save(fullfile(folderImgTrData,'DataTrPatchFeat.mat'),'trFeatMat');
else
    load(fullfile(folderImgTrData,'DataTrPatchFeat.mat'),'trFeatMat');
end

%% Train a patch (bolt) classifier
if ~exist(fullfile(folderImgTrData,'DataTrClss.mat'),'file');
    
    % Gentle boost training.
    clss = fitensemble(trFeatMat,[boltGTPatch(:).clss],'GentleBoost', ...
        nRounds, 'Tree','Type','Classification');

    weakIdx = zeros(nRounds,1);
    for ii = 1:nRounds
        % There are "x" in the begging of index names. 
        tmp = clss.Trained{ii}.CompactRegressionLearner.CutPredictor{1,1};
        weakIdx(ii) = str2double(tmp(2:end));
    end
    
    % useful features from all candidate features
    weakIdx = unique(weakIdx); 
   
    % This process complecate is that we extract a (nHaar*nChannel) number
    % of features from images however, they don't contribute the decision
    % making (weak classifier). Then, we don't need to extract all
    % nHaar*nChannel features.
    
    % Re-training for reducing 'clss' dimension using a valid features,
    % which are used for the above training.
    clss = fitensemble(trFeatMat(:,weakIdx),[boltGTPatch(:).clss],...
        'GentleBoost', nRounds, 'Tree','Type','Classification');
    
    % compact the trained classifier
    clss = compact(clss);
    save(fullfile(folderImgTrData,'DataTrClss'),'clss','weakIdx');
else
    load(fullfile(folderImgTrData,'DataTrClss'),'clss','weakIdx'); 
end

%% Hard negative training: patch detection on training images
if ~exist(fullfile(folderImgTrData,'DataTrPatchDetect.mat'),'file');
    nTrImg      = numel(rszTrImgs);
    detectPatch = cell(nTrImg,nScale);
    for ii=1:nTrImg
        img     = imread(rszTrImgs{ii}); % load a training image
        for jj  = 1: nScale
            imgR        = imresize(img, imgScale(jj));
            imsz        = [size(imgR,1) size(imgR,2)];
            
            itgImgSet   = CompIntChnImg(imgR);
            
            % This is very useful function of producing overlap sliding
            % windows found in https://github.com/tsogkas/matlab-utils
            % Users can use "blockproc" in Matlab but, I don't know how to
            % make this work with having stride. For example, moving the
            % sliding windows with 4 pixel strides. 
            [~,ind,~]   =  im2patches( ...
                reshape(1:(imsz(1)+1)*(imsz(2)+1),imsz(1)+1, ...
                imsz(2)+1), [patchSize+1,patchSize+1], patchStride);
            
            % Let's only consider fully included patches. Remove patches
            % overlapping with image boundary.
            centIdx = sub2ind([patchSize+1 patchSize+1], ...
                patchSize/2+1, patchSize/2+1);
            
            [y,x] = ind2sub([imsz(1)+1, imsz(1)+1], ind(centIdx, :));
            
            validIdx = (y >= patchSize/2+1) & ...
                       (y <= imsz(1) - (patchSize/2+1)) & ...
                       (x >= patchSize/2+1) & ...
                       (x <= imsz(2) - (patchSize/2+1));
            
            ind = ind(:,validIdx);

%             % slow version
%             for kk=1:nWeakClss
%                 haarComp = haarPatch{featIdx.haar(weakIdx(kk)),2};
%                 itgImg   = itgImgSet(:,:,featIdx.channel(weakIdx(kk)));
%                 [~,imcol,~] =  im2patches(itgImg, ...
%                     [patchSize+1,patchSize+1], patchStride);
%                 featImgSet(:,kk) = haarComp(:)'*imcol;  
%             end; clearvars imcol;

            % the above "slow" version is easy to understand how it works.
            % faster version (exactly process)
            featImgSet  = zeros(size(ind,2),numel(weakIdx));
            
            haarList = featIdx.haar(weakIdx(:));
            chanList = featIdx.channel(weakIdx(:));
            for kk=1:nChannel
                itgImg   = itgImgSet(:,:,kk);
                [~,imcol,~] =  im2patches(itgImg, ...
                    [patchSize+1,patchSize+1], patchStride);
                
                imcol = imcol(:,validIdx);
                
                haarChan = find(chanList==kk);
                for qq=1:numel(haarChan)
                    haarComp = haarPatch{haarList(haarChan(qq)),2};
                    featImgSet(:,haarChan(qq)) = haarComp(:)'*imcol;
                end
            end

            clssOut = predict(clss,featImgSet);
            
            % the size of the detected patch is same.
            % x and y below are the center of each detected windows. 
            % positive samples
            [y,x] = ind2sub([imsz(1)+1, imsz(1)+1], ...
                ind(centIdx, find(clssOut==1)));
            nDetect = numel(x);
            detectPatch{ii,jj,1} = [x',y', ...
                ones(nDetect,1)*patchSize, ...
                ones(nDetect,1)*patchSize]./imgScale(jj);
            
            % random negative samples: nNegMul is to generate similar
            % number of random negative samples with the positive samples
            % for new training. 
            [y,x] = ind2sub([imsz(1)+1, imsz(1)+1], ...
                ind(centIdx, find(clssOut==-1)));
            sampleNeg  = randperm(numel(x),nDetect*nNegMul);
            
            detectPatch{ii,jj,2} = [x(sampleNeg)',y(sampleNeg)', ...
                ones(nDetect*nNegMul,1)*patchSize, ...
                ones(nDetect*nNegMul,1)*patchSize]./imgScale(jj);
            
            fprintf('Processing %d(/%d) image of %d(/%d) scale \n', ...
                ii, nTrImg,jj,nScale);
        end
    end
    save(fullfile(folderImgTrData,'DataTrPatchDetect'),'detectPatch');
else
    load(fullfile(folderImgTrData,'DataTrPatchDetect'),'detectPatch');
end

%% Hard negative training: extraction of hard negative Patch 
if ~exist(fullfile(folderImgTrData,'DataTrValidPatch.mat'),'file');
    nTrImg = numel(rszTrImgs);
    boltNewPatch = cell(nTrImg,2);
    for ii=1:nTrImg
        
        detectPosRect = [];
        detectNegRect = [];
        for jj  = 1: nScale
            detectPosRect = [detectPosRect;detectPatch{ii,jj,1}];
            detectNegRect = [detectNegRect;detectPatch{ii,jj,2}];
        end
        
        % "rect" order
        detectPosRect(:,1:2) = detectPosRect(:,1:2)-detectPosRect(:,3:4)/2;
        detectNegRect(:,1:2) = detectNegRect(:,1:2)-detectNegRect(:,3:4)/2;
         
        idxPos   = find(and([boltGTPatch(:).trImgIdx]==ii, ...
            [boltGTPatch(:).clss] == 1));
        idxNeg   = find(and([boltGTPatch(:).trImgIdx]==ii, ...
            [boltGTPatch(:).clss] == -1));
        
        % True-psotive and false-negative samples in the previous step are
        % included in positive and negative samples for new training,
        % repectively.
        
        % ground-truth positive boxes
        bboxA = reshape([boltGTPatch(idxPos).rszOrgBbox],4,numel(idxPos));
        bboxA = permute(bboxA,[2 1]);
        
        % positive detection
        bboxB = detectPosRect;
        
        % IOU computation between positive detection and GT boxes.
        tmp = bboxOverlapRatio(bboxA,bboxB);
        boltNewPatch{ii,1} = ...
            [detectPosRect(logical(sum(tmp)>thresh_IOU),:);bboxA];
        
        % ground-truth negative boxes
        bbox = reshape([boltGTPatch(idxNeg).rszOrgBbox],4,numel(idxNeg));
        bbox = permute(bbox,[2 1]);
        boltNewPatch{ii,2} = ...
        [detectPosRect(logical(sum(tmp)<thresh_IOU),:);bbox;detectNegRect];
    end
    save(fullfile(folderImgTrData,'DataTrValidPatch'),'boltNewPatch');
else
    load(fullfile(folderImgTrData,'DataTrValidPatch'),'boltNewPatch');
end

%% Train a patch classifier
if ~exist(fullfile(folderImgTrData,'DataTrNewPatchFeat.mat'),'file');

    nTrImg = numel(rszTrImgs);
    clssPatch       = cell(nTrImg,1);
    trNewFeatMat    = cell(nTrImg,1);
    for ii=1:nTrImg
        % load a training image
        img     = imread(rszTrImgs{ii});
        
        nPos    = size(boltNewPatch{ii,1},1);
        nNeg    = size(boltNewPatch{ii,2},1);
        clssPatch{ii} = [ones(nPos,1);ones(nNeg,1)*-1];
        
        % positive
        for jj=1:nPos
            boltPatchTmp = ...
                imresize(imcrop(img,boltNewPatch{ii,1}(jj,:)), ...
                [patchSize patchSize]);
            
            % compute features from integral channel images
            trNewFeatMat{ii}(jj,:) = ...
                CompFeatChnPatch(CompIntChnImg(boltPatchTmp), ...
                haarPatch, nHaar);
        end
        
        % negative
        for jj=1:nNeg
            boltPatchTmp = ...
                imresize(imcrop(img,boltNewPatch{ii,2}(jj,:)), ...
                [patchSize patchSize]);
            
            % compute features from integral channel images
            trNewFeatMat{ii}(jj+nPos,:) = ...
                CompFeatChnPatch(CompIntChnImg(boltPatchTmp), ...
                haarPatch, nHaar);
        end
    end
    
    trNewFeatMatAll  = [];
    trClssPatch = [];
    for ii=1:nTrImg
        trNewFeatMatAll = [trNewFeatMatAll;trNewFeatMat{ii}];
        trClssPatch     = [trClssPatch;clssPatch{ii}];
    end; clearvars trNewFeatMat clssPatch;
    
    clss = fitensemble(trNewFeatMatAll,trClssPatch,'GentleBoost', ...
        nRounds, 'Tree','Type','Classification');
    
    weakNewIdx = zeros(nRounds,1);
    for ii = 1:nRounds
        tmp = clss.Trained{ii}.CompactRegressionLearner.CutPredictor{1,1};
        weakNewIdx(ii) = str2double(tmp(2:end));
    end
    weakNewIdx = unique(weakNewIdx);
    
    % re-training for reducing 'clss' dimension
    clss = fitensemble(trNewFeatMatAll(:,weakNewIdx), ...
        trClssPatch,'GentleBoost', nRounds, ...
        'Tree','Type','Classification');
    
    clss = compact(clss);
    save(fullfile(folderImgTrData,'DataTrFinalClss'),'clss','weakNewIdx');
else
    load(fullfile(folderImgTrData,'DataTrFinalClss'),'clss','weakNewIdx');
end