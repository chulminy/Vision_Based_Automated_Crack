Parameters;
%% Step 2: Object detect
% load training data
if ~exist(fullfile(folderImgTsData,'DataImg.mat'),'file');
    
    load(fullfile(folderImgTrData,'DataTrFinalClss.mat'), ...
        'clss','weakNewIdx');
    load(fullfile(folderImgTrData,'DataHaar.mat'),'haarPatch','featIdx');

    nTsImg  = numel(rszTsImgs);
    
    % initializing DataImg
    for ii=1:nTsImg
        DataImg(ii).filename = rszTsImgs{ii};
    end
    
    for ii=1:nTsImg
        % image resizing --------------------------------------------------
        if ~exist(rszTsImgs{ii},'file')
            img = imresize(imread(orgTsImgs{ii}),imgRszScale);
            img = imfilter(img, gaussF); % remove noise
            imwrite(img,rszTsImgs{ii});
        else
            img = imread(rszTsImgs{ii});
        end
                
        % object detection ------------------------------------------------
        detectRect = cell(nScale,1);
        for jj  = 1: nScale
            imgR        = imresize(img, imgScale(jj));
            imsz        = [size(imgR,1) size(imgR,2)];
            
            itgImgSet   = CompIntChnImg(imgR);
            
            [~,ind,~]   = im2patches( ...
                reshape(1:size(itgImgSet,1)*size(itgImgSet,2), ...
                size(itgImgSet,1), size(itgImgSet,2)),...
                [patchSize+1,patchSize+1], patchStride);
            
            centIdx = sub2ind([patchSize+1 patchSize+1], ...
                patchSize/2+1, patchSize/2+1);
            
            [y,x] = ind2sub([imsz(1)+1, imsz(1)+1], ind(centIdx, :));
            
            validIdx = (y >= patchSize/2+1) & ...
                       (y <= imsz(1) - (patchSize/2+1)) & ...
                       (x >= patchSize/2+1) & ...
                       (x <= imsz(2) - (patchSize/2+1));
            
            ind = ind(:,validIdx);
            
            featImgSet  = zeros(size(ind,2),numel(weakNewIdx));
            
            haarList = featIdx.haar(weakNewIdx(:));
            chanList = featIdx.channel(weakNewIdx(:));
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
            
            [y,x]   = ind2sub([size(itgImgSet), size(itgImgSet,2)], ...
                ind(centIdx, find(clssOut==1)));
            
            % note that x and y are the center of each patch.
            
            detectRect{jj} = [x',y', ...
                ones(numel(x),1)*patchSize, ...
                ones(numel(x),1)*patchSize]./imgScale(jj);
           
            fprintf('Processing %d(/%d) test image of %d(/%d) scale \n',...
                ii,nTsImg,jj,nScale);
        end
        
        detectBolt = [];
        for jj  = 1: nScale
            detectBolt = [detectBolt;(detectRect{jj})];
        end
        
        % "rect" order
        detectBolt(:,1:2) = detectBolt(:,1:2)-detectBolt(:,3:4)/2;
        
        % merging the detected windows
        overlapRatio  = bboxOverlapRatio(detectBolt,detectBolt);
        
        % construct adj matrix
        overlapRatioAdj = overlapRatio>thresh_IOU;
        
        % apply community detection
        com = GCModulMax3(overlapRatioAdj);
        
        % merging
        comIdx = unique(com);
        objFinal = zeros(numel(comIdx),4);
        for jj=1:numel(comIdx)
            tmp = detectBolt(find(com==jj),:);
            % the bounding box encompassing all valid patches
            tmp(:,3:4) = tmp(:,1:2) + tmp(:,3:4);
            objFinal(jj,1:2) = min(tmp(:,1:2),[],1);
            objFinal(jj,3:4) = max(tmp(:,3:4),[],1); 
        end
        objFinal(:,3:4) = objFinal(:,3:4)-objFinal(:,1:2);
        
        nObj      = size(objFinal,1);
        
        featurePt = DataImg(ii).F(1:2,:);
        DataImg(ii).nObj        = nObj;
        for jj=1:nObj
            points  = bbox2points(objFinal(jj,:));
            in      = inpolygon(featurePt(1,:)',featurePt(2,:)', ...
                points(:,1),points(:,2));
            
            DataImg(ii).obj(jj).bbox = objFinal(jj,:);
            DataImg(ii).obj(jj).cPt  = objFinal(jj,1:2)+objFinal(jj,3:4)/2;
            DataImg(ii).obj(jj).featInObj =  in;
        end
        
        % sift feature extraction for object grouping ---------------------
        [F, D] = vl_sift(single(rgb2gray(img)));
        DataImg(ii).F  = F;
        DataImg(ii).D  = D;
    end
    save(fullfile(folderImgTsData,'DataImg.mat'),'DataImg');
else
    load(fullfile(folderImgTsData,'DataImg.mat'),'DataImg');
end

%% Step 3:  Object grouping
% compute fundamental matrix ----------------------------------------------
if ~exist(fullfile(folderImgTsData,'DataFunMat.mat'),'file');
    % assumption that images are captured consecutively.
    % fundamental matrixes are computed next numCompImg images.
    
    % find SIFT features
    nTsImg = numel(rszTsImgs);
 
    fundMat     = cell(nTsImg,nTsImg);
    for ii=1:nTsImg
        F1 = DataImg(ii).F;   
        D1 = DataImg(ii).D;      
        for jj=ii+1:min(numCompImg+ii,nTsImg)
            F2 = DataImg(jj).F;   
            D2 = DataImg(jj).D;    
            
            % VL_FEAT matching
            [matches, ~] = vl_ubcmatch(D1,D2);
            
            X1 = F1(1:2,matches(1,:));
            X2 = F2(1:2,matches(2,:));
            
            [fund,inliers,status] = estimateFundamentalMatrix(X1',X2');
            
            if status ~= 0;
                fundMat{ii,jj} = eye(3,3);
            else
                fundMat{ii,jj}.fund   = fund;
                fundMat{ii,jj}.matches = matches(:,inliers);
            end;
        end
    end
    save(fullfile(folderImgTsData,'DataFunMat.mat'),'fundMat');
else
    load(fullfile(folderImgTsData,'DataFunMat.mat'),'fundMat');
end

% generate DataObj --------------------------------------------------------
if ~exist(fullfile(folderImgTsData,'DataObj.mat'),'file');
    nTsImg  = numel(rszTsImgs);
    count   = 1;
    for ii=1:nTsImg
        nObj = DataImg(ii).nObj;
        for jj=1:nObj
            DataObj(count).imgIdx    = ii;
            DataObj(count).bbox      = DataImg(ii).obj(jj).bbox;
            DataObj(count).cPt       = DataImg(ii).obj(jj).cPt;
            DataObj(count).featInObj = DataImg(ii).obj(jj).featInObj;
            count = count + 1;
        end
    end
    save(fullfile(folderImgTsData,'DataObj.mat'),'DataObj');
else
    load(fullfile(folderImgTsData,'DataObj.mat'),'DataObj');
end

if ~exist(fullfile(folderImgTsData,'DataBolt.mat'),'file');
    nTotalObj = numel(DataObj);
    
    adjMat = false(nTotalObj,nTotalObj);
    matchCountMat = zeros(nTotalObj,nTotalObj);
    for ii=1:nTotalObj
        img1Idx = DataObj(ii).imgIdx;
        for jj=ii+1:nTotalObj
            
            img2Idx = DataObj(jj).imgIdx;
            if img2Idx == img1Idx; continue; end;
            if img2Idx > img1Idx + numCompImg; continue; end;
            
            F = fundMat{img1Idx,img2Idx}.fund;
            
            matches = fundMat{img1Idx,img2Idx}.matches;
            
            [~,Lb1] = ismember(find(DataObj(ii).featInObj), matches(1,:));
            [~,Lb2] = ismember(find(DataObj(jj).featInObj), matches(2,:));
            
            featMatch = nonzeros(intersect(Lb1,Lb2));
            
            if numel(featMatch) > nMinFeatMatch;  
                adjMat(ii,jj) = true; 
                matchCountMat(ii,jj) = numel(featMatch);
            end

        end
        % apply a constraint that a bolt in a image must be corresponded to one
        % in another images. "matchCountMat" is used.
        adjTs    = adjMat(ii,:);
        countMat = matchCountMat(ii,:);
        
        imgIdx   = [DataObj(:).imgIdx];
        for jj=1:nTsImg
            idxTs = find(imgIdx==jj);
            if sum(adjTs(idxTs))> 1
                [~,id] = max(countMat(idxTs));
                adjMat(ii,idxTs) = false;
                adjMat(ii,idxTs(id)) = true;
            end
        end
    end
    
    % apply community detection
    com = GCModulMax3(adjMat);
    
    % remove communities having small number of memebers.
    comIdx = unique(com);
    objList = [];
    for jj=1:numel(comIdx)
        objIdx = find(com==jj);
        nPatch = numel(objIdx);
        if nPatch < nMinPatch; 
            objList = [objList;objIdx];
        end   
    end
    adjMat(objList,:) = [];
    adjMat(:,objList) = [];
    DataObj(objList)  = [];
    
    % re-apply community detection
    com = GCModulMax3(adjMat);
    
    % merging
    comIdx = unique(com);
    imgIdx = [DataObj(:).imgIdx];
    bbox   = cell2mat({DataObj(:).bbox}');
    
    count  = 1;
    for jj=1:numel(comIdx)
        nPatch = numel(imgIdx(find(com==jj)));
        if nPatch < nMinPatch; continue; end;
        DataBolt(count).imgIdx = imgIdx(find(com==jj))';
        DataBolt(count).bbox   = bbox(find(com==jj),:);
        DataBolt(count).nPatch = nPatch;
        count = count + 1;
    end
    save(fullfile(folderImgTsData,'DataBolt.mat'),'DataBolt');
else
    load(fullfile(folderImgTsData,'DataBolt.mat'),'DataBolt');
end

%% Step 4:  Crack Detection
if ~exist(fullfile(folderImgTsData,'DataCrack.mat'),'file');
    
nBolt = numel(DataBolt);
count = 1;
for ii=1:nBolt

    imgIdx      = DataBolt(ii).imgIdx;
    bbox        = DataBolt(ii).bbox;
    nBoltImg    = DataBolt(ii).nPatch;
    for jj=1:nBoltImg
        
        % image read
        img     = imread(orgTsImgs{imgIdx(jj)});
        imgSizeOrg = [size(img,1) size(img,2)];
        
        % regularizing object patch
        bboxBolt = bbox(jj,:)./imgRszScale;
        
        cent    = bboxBolt(1:2)+bboxBolt(3:4)/2;
        winS    = max(bboxBolt(3:4))*scaleMul;
        
        if (cent(1)-winS/2) < 1 || (cent(2)-winS/2 < 1)
            winS = 2 * min(cent(1)-1, cent(2)-1);
        end
        
        if (cent(1)+winS/2) > imgSizeOrg(2) || ...
                (cent(2)+winS/2 > imgSizeOrg(1))
            winS = 2 * min(imgSizeOrg(2)-cent(1), imgSizeOrg(1)-cent(2));
        end
        
        bboxBolt = [fix(cent - [winS/2 winS/2]) winS winS];
        
        
        patch   = im2double(rgb2gray(imresize(imcrop(img,bboxBolt),...
            [patchSizeOrg * scaleMul patchSizeOrg * scaleMul])));
        
        [crackBW, flag] = CrackDetection(patch, Param);
        
        if flag == 1;
            orgBW = imresize(crackBW,bboxBolt(3:4));
            stats = regionprops(orgBW,'PixelList');
            
            if numel(stats)>1;
                DataCrack(count).boltID = ii;
                DataCrack(count).imgIdx = imgIdx(jj);
                DataCrack(count).bbox   = bboxBolt;
                tmp = [];
                for qq=1:numel(stats)
                 tmp = [tmp; ...
                    bsxfun(@plus,[stats(qq).PixelList],bboxBolt(1:2))];
                end
                DataCrack(count).pxList = tmp
                count = count + 1;
                
            else
                DataCrack(count).boltID = ii;
                DataCrack(count).imgIdx = imgIdx(jj);
                DataCrack(count).bbox   = bboxBolt;
                DataCrack(count).pxList = ...
                    bsxfun(@plus,[stats.PixelList],bboxBolt(1:2));
                count = count + 1;
            end
        end
    end
end
    save(fullfile(folderImgTsData,'DataCrack.mat'),'DataCrack');
else
    load(fullfile(folderImgTsData,'DataCrack.mat'),'DataCrack');
end


