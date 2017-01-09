% This file is to plot intermediate and final outcomes. You should first
% run RunTrainClassifier.m and RunTest.m. You need to activate the
% following flags to selectively plot the outcome(s). All plotted images
% will be saved in folderOut.
Parameters;

% Set "Yes" in the outcome that you want.

% (training) Show original training images
PLOTOUT1 = false;

% (training) Show positive and negative patches for training
PLOTOUT2 = false;

% (training) Show channel patches
PLOTOUT3 = false;

% (training) Show classified patches on training images
PLOTOUT4 = false;

% (testing) Show testing images
PLOTOUT5 = false;

% (testing) Show detected objects
PLOTOUT6 = false;

% (testing) Show a group of each bolt
PLOTOUT7 = false;

% (testing) Show crack detection results
PLOTOUT8 = false;

%% PLOTOUT1: (training) Show training images
if PLOTOUT1
    nTrImg  = numel(rszTrImgs); % # of training images
    imgSize = size(imresize(imread(rszTrImgs{1}),0.3));
    imageset = zeros(imgSize(1),imgSize(2),3,nTrImg,'uint8');
    prm = struct('mm',2,'nn',3,'padAmt',1,'hasChn',1,'showLines',1);
    for ii=1:nTrImg
        imageset(:,:,:,ii) = imresize(imread(rszTrImgs{ii}),0.3);
    end
    figure('Name','PlotOut1'); clf; h = montage2(imageset, prm );
    imwrite(h.CData,fullfile(folderOut,'PlotOut1.jpg'));
    imwrite(imresize(h.CData,[NaN 900]), ...
        fullfile(folderHome,'PlotOut1.jpg'));
    clearvars imageset imgSize nTrImg prm
end

%% PLOTOUT2: (training) Show training images
if PLOTOUT2
    load(fullfile(folderImgTrData,'DataTrGTPatch.mat'),'boltGTPatch');
    prm = struct('padAmt',1,'hasChn',1,'showLines',1);
    posIdx = find([boltGTPatch(:).clss]==1);
    tmp = {boltGTPatch(posIdx(randperm(numel(posIdx),50))).patch};
    figure(1);subplot(121);h1=montage2(cat(4,tmp{:}),prm); clearvars tmp;
    negIdx = find([boltGTPatch(:).clss]==-1);
    tmp = {boltGTPatch(negIdx(randperm(numel(negIdx),50))).patch};
    subplot(122);h2=montage2(cat(4,tmp{:}),prm); clearvars tmp prm;
    
    img = cat(2,h1.CData,h2.CData);
    imwrite(img,fullfile(folderOut,'PlotOut2.jpg'));
    imwrite(imresize(img,[NaN 900]),fullfile(folderHome,'PlotOut2_2.jpg'));
    clearvars img prm
end

%% PLOTOUT3: (training) Show channel patches
if PLOTOUT3
    load(fullfile(folderImgTrData,'DataTrGTPatch.mat'),'boltGTPatch');
    idBoltGTPatch = 5;
    
    img = boltGTPatch(idBoltGTPatch).patch;
    % LUV transformation---------------------------------------------------
    patchLUV = rgbConvert(img, 'luv');
    patchHSV = rgbConvert(img, 'hsv');
    
    % Gradient histogram---------------------------------------------------
    [M,O]       = gradientMag(single(img));
    nOrients    = 6;
    H           = gradientHist(M,O,1,nOrients,0,0,0);
    
    imgSet1 = mat2gray(rgb2gray(img));
    imgSet1 = cat(2,imgSet1,mat2gray(patchLUV(:,:,2)));
    imgSet1 = cat(2,imgSet1,mat2gray(patchLUV(:,:,3)));
    imgSet1 = cat(2,imgSet1,mat2gray(patchHSV(:,:,2)));
    imgSet1 = cat(2,imgSet1,mat2gray(patchHSV(:,:,3)));
    imgSet1 = cat(2,imgSet1,mat2gray(M));
    
    imgSet2 = mat2gray(H(:,:,1));
    imgSet2 = cat(2,imgSet2,mat2gray(H(:,:,2)));
    imgSet2 = cat(2,imgSet2,mat2gray(H(:,:,3)));
    imgSet2 = cat(2,imgSet2,mat2gray(H(:,:,4)));
    imgSet2 = cat(2,imgSet2,mat2gray(H(:,:,5)));
    imgSet2 = cat(2,imgSet2,mat2gray(H(:,:,6)));
    
    imgSet  = cat(1,imgSet1,imgSet2);
    
    figure('Name','PlotOut1'); imshow(imgSet);
    imwrite(imgSet,fullfile(folderOut,'PlotOut3.jpg'));
    imwrite(imresize(imgSet,[NaN 900]), ...
        fullfile(folderHome,'PlotOut3.jpg'));
    clearvars imgSet imgSet1 imgSet2 H M O
end

%% PLOTOUT4: (training) Show classified patches on training images
if PLOTOUT4
    load(fullfile(folderImgTrData,'DataTrPatchDetect'),'detectPatch');
    load(fullfile(folderImgTrData,'DataTrGTPatch.mat'),'boltGTPatch');
    indTrImg = 1;
    
    indx  = intersect(find([boltGTPatch(:).trImgIdx]==indTrImg), ...
        find([boltGTPatch(:).clss]==1));
    
    GTBox = reshape(cell2mat({boltGTPatch(indx).rszOrgBbox}),4,numel(indx));
    detectPosRect = [];
    detectNegRect = [];
    for jj  = 1: nScale
        detectPosRect = [detectPosRect;detectPatch{indTrImg,jj,1}];
        detectNegRect = [detectNegRect;detectPatch{indTrImg,jj,2}];
    end
    
    % "rect" order
    detectPosRect(:,1:2) = detectPosRect(:,1:2)-detectPosRect(:,3:4)/2;
    detectNegRect(:,1:2) = detectNegRect(:,1:2)-detectNegRect(:,3:4)/2;
    
    nPos = size(detectPosRect,1);
    nNeg = size(detectNegRect,1);
    
    img  = imread(rszTrImgs{indTrImg});
    
    img = insertShape(img,'Rectangle', ...s
        detectNegRect(randperm(nNeg,50),:), ...
        'LineWidth',10,'Color','red');
    
    img = insertShape(img,'Rectangle', ...
        detectPosRect(randperm(nPos,50),:), ...
        'LineWidth',10,'Color','blue');
    
    img = insertShape(img,'Rectangle', GTBox', ...
        'LineWidth',10,'Color','green');
    
    figure('Name','PlotOut4'); imshow(img);
    imwrite(img,fullfile(folderOut,'PlotOut4.jpg'));
    imwrite(imresize(img,[NaN 900]), ...
        fullfile(folderHome,'PlotOut4.jpg'));
end

%% PLOTOUT5: (testing) Show testing images
if PLOTOUT5
    nTsImg  = numel(rszTsImgs); % # of training images
    imgSize = size(imresize(imread(rszTsImgs{1}),0.3));
    imageset = zeros(imgSize(1),imgSize(2),3,nTsImg,'uint8');
    prm = struct('mm',5,'nn',14,'padAmt',1,'hasChn',1,'showLines',1);
    for ii=1:nTsImg
        imageset(:,:,:,ii) = imresize(imread(rszTsImgs{ii}),0.3);
    end
    figure('Name','PlotOut1'); clf; h = montage2(imageset, prm );
    imwrite(h.CData,fullfile(folderOut,'PlotOut5.jpg'));
    imwrite(imresize(h.CData,[NaN 900]), ...
        fullfile(folderHome,'PlotOut5.jpg'));
    clearvars imageset imgSize nTsImg prm
end

%% PLOTOUT6: (testing) Show detected objects
if PLOTOUT6
    load(fullfile(folderImgTsData,'DataImg.mat'),'DataImg');
    
    indTsImg = randperm(numel(rszTsImgs),6);
    imgSet   = cell(6,1);
    for ii=1:numel(indTsImg)
        img  = imread(rszTsImgs{indTsImg(ii)});
        nObj = DataImg(indTsImg(ii)).nObj;
        detectedBox = [];
        for jj=1:nObj
            detectedBox = [detectedBox;DataImg(indTsImg(ii)).obj(jj).bbox];
        end
        imgSet{ii} = insertShape(img,'Rectangle', detectedBox, ...
            'LineWidth',10,'Color','blue');
    end
    
    img1 = imgSet{1};
    img1 = cat(2,img1,imgSet{2});
    img1 = cat(2,img1,imgSet{3});
    
    img2 = imgSet{4};
    img2 = cat(2,img2,imgSet{5});
    img2 = cat(2,img2,imgSet{6});
    
    img  = cat(1,img1,img2);
    
    figure('Name','PLOTOUT6'); imshow(img);
    imwrite(img,fullfile(folderOut,'PLOTOUT6.jpg'));
    imwrite(imresize(img,[NaN 900]), ...
        fullfile(folderHome,'PLOTOUT6.jpg'));
    clearvars imageset imgSize nTsImg prm
end

%% PLOTOUT7: (testing) Show a group of each bolt
if PLOTOUT7
    load(fullfile(folderImgTsData,'DataBolt.mat'),'DataBolt');
    
    nBolt = numel(DataBolt);
    
    nPatch  = 0;
    while nPatch < 16
        idxBolt = randperm(nBolt,1);
        nPatch = DataBolt(idxBolt).nPatch;
    end;
    imgIdx = DataBolt(idxBolt).imgIdx;
    bbox   = DataBolt(idxBolt).bbox;
    
    imgSet   = cell(16,1);
    for ii=1:16
        img  = imread(rszTsImgs{imgIdx(ii)});
        img  = insertShape(img,'Rectangle', bbox(ii,:), ...
            'LineWidth',10,'Color','blue');
        imgSet{ii} = imresize(img,0.3);
    end
    
    patchSet = cell(16,1);
    for ii=1:16
        img  = imread(rszTsImgs{imgIdx(ii)});
        patch   = imcrop(img,bbox(ii,:));
        patchSet{ii} = patch;
    end
    
    imgAll   = [];
    patchAll = [];
    count = 1;
    for ii=1:4
        imgTmp = [];
        patchTmp = [];
        for jj=1:4
            imgTmp = cat(2,imgTmp,imgSet{count});
            patchTmp = cat(2,patchTmp, ...
                imresize(patchSet{count},[patchSizeOrg patchSizeOrg]));
            count = count + 1;
        end
        imgAll      = cat(1,imgAll,imgTmp);
        patchAll    = cat(1,patchAll,patchTmp);
    end
    
    % insert space
    imgTmp = cat(2,imresize(patchAll,[size(imgAll,1) NaN]), ...
        imresize(ones(1,1,3)*255,[size(imgAll,1) 300]));
    
    img = cat(2,imgTmp, imgAll);
    
    figure('Name','PLOTOUT7'); imshow(img);
    imwrite(img,fullfile(folderOut,'PLOTOUT7.jpg'));
    imwrite(imresize(img,[NaN 900]), ...
        fullfile(folderHome,'PLOTOUT7.jpg'));
    clearvars imageset imgSize nTsImg prm
end

%% % (testing) Show crack detection results
if PLOTOUT8
    load(fullfile(folderImgTsData,'DataCrack.mat'),'DataCrack');
    
    boltID = [DataCrack(:).boltID];
    
    crackedBolt  = unique(boltID);
    nCrackedBolt = numel(crackedBolt);
    
    imgAll = [];
    
    for ii=1:nCrackedBolt
        idxPatch = find(boltID==crackedBolt(ii));
        
        imgAll = [];
        for jj=1:numel(idxPatch)
            
            imgTmp = [];
            imgIdx = DataCrack(idxPatch(jj)).imgIdx;
            bbox   = DataCrack(idxPatch(jj)).bbox;
            
            img    = imread(orgTsImgs{imgIdx});
            patch1 = imcrop(img,bbox);
            
            imgC   = img;
            pxList = DataCrack(idxPatch(jj)).pxList;
            nPx = size(pxList,1);
            for kk=1:nPx
                imgC(pxList(kk,2),pxList(kk,1),:) = [255 0 0];
            end
            patch2 = imcrop(imgC,bbox);
            
            img  = insertShape(img,'Rectangle', bbox, ...
                'LineWidth',20,'Color','blue');
            
            imgC  = insertShape(imgC,'Rectangle', bbox, ...
                'LineWidth',20,'Color','blue');
            
            img1   = imresize(img,0.1);
            img2   = imresize(imgC,0.1);
            imSize = [size(img1,1) size(img1,2)];
            
            imgTmp = cat(2,imgTmp,imresize(patch1,[imSize(1) NaN]));
            imgTmp = cat(2,imgTmp,imresize(patch2,[imSize(1) NaN]));
            imgTmp = cat(2,imgTmp,img1);
            imgTmp = cat(2,imgTmp,img2);
            
            [~,name,ext] = fileparts(orgTsImgs{imgIdx});
            
            imgTmp = insertText(imgTmp,[ceil(size(imgTmp,2)/2)+10 10], ...
                ['Bolt ID: ' int2str(crackedBolt(ii)) ...
                ' and Image ID: ' name ext], ...
                'FontSize',25, 'Font','Arial Bold');
            
            if jj == numel(idxPatch);
                imgAll = cat(1,imgAll,imgTmp);
            else
                imgAll = cat(1,imgAll,imgTmp);
                imgAll = cat(1,imgAll,ones(5,size(imgTmp,2),3)*255);
            end
        end
        
        figure('Name',['PLOTOUT8-BoltID: ' int2str(crackedBolt(ii))]); 
        imshow(imgAll); imwrite(imgAll,fullfile(folderOut, ...
            ['PLOTOUT8-BoltID-' int2str(crackedBolt(ii)) '.jpg']));
        imwrite(imresize(imgAll,[NaN 900]), ...
            fullfile(folderHome, ...
            ['PLOTOUT8-BoltID-' int2str(crackedBolt(ii)) '.jpg']));
        
    end
end

% -------------------------------------------------------------------------
% Below code is to create images used for web posting. Just ignore.
webPost = true;
if webPost
    Parameters;
    
    load(fullfile(folderImgTsData,'DataCrack.mat'),'DataCrack');
    
    boltID = [DataCrack(:).boltID];
    
    crackedBolt  = unique(boltID);
    nCrackedBolt = numel(crackedBolt);
    
    imgAll = [];
    
    for ii=5:5
        idxPatch = find(boltID==crackedBolt(ii));
        
        imgAll = [];
        for jj=2:2
            
            imgTmp = [];
            imgIdx = DataCrack(idxPatch(jj)).imgIdx;
            bbox   = DataCrack(idxPatch(jj)).bbox;
            
            img    = imread(orgTsImgs{imgIdx});
            patch1 = imcrop(img,bbox);
            
            imgC   = img;
            pxList = DataCrack(idxPatch(jj)).pxList;
            nPx = size(pxList,1);
            for kk=1:nPx
                imgC(pxList(kk,2),pxList(kk,1),:) = [255 0 0];
            end
            patch2 = imcrop(imgC,bbox);
            
            img  = insertShape(img,'Rectangle', bbox, ...
                'LineWidth',20,'Color','blue');
            
            imgC  = insertShape(imgC,'Rectangle', bbox, ...
                'LineWidth',20,'Color','blue');
            
            img1   = imresize(img,0.1);
            img2   = imresize(imgC,0.1);
            imSize = [size(img1,1) size(img1,2)];
            
            imgTmp = []; 
            imgTmp = cat(2,imgTmp,imresize(patch1,[imSize(1) NaN]));
            imgTmp = cat(2,imgTmp,imresize(patch2,[imSize(1) NaN]));
           
            [~,name,ext] = fileparts(orgTsImgs{imgIdx});
            
        end

        imwrite(imgTmp,fullfile(folderHome,'teaser.jpg'));
    end
end


