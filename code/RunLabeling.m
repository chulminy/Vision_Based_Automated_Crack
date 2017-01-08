%% Description 
%
% This tool is to annotate bounding boxes (patches) which include an
% objects of interest, a bolts here. All images in the training image
% folder should be labaeld. Once the annotation files (["DataLabel" idx
% .m]) of all training images are generated, users need to run the
% following code to merge the annotation data: Let's say the number of
% training images are "nImg" here.
%
%{
Parameters;
for ii=1:5
    load(fullfile(folderImgTrData,['DataLabel' int2str(ii)]));
    boltTrLabel(ii).filename   = orgTrImgs{ii};
    boltTrLabel(ii).bbox       = bbox;
    boltTrLabel(ii).clssLabel  = classLabel;
end
save(fullfile(folderImgTrData,'DataLabelAll'),'boltTrLabel');
%}
Parameters;

% Users need to update this!!
%--------------------------------------------------------------------------
idx = input('Please type the image index: '); 
%--------------------------------------------------------------------------

% show the image
fig = figure(1); 
img = imread(orgTrImgs{idx}); imshow(img);
title(['Draw a rectangular with a mouse.' ...
       'Press any to exit the selection'],'FontSize',14);axis equal;

% mask image (positive region = 1, other. 0) to find negative patches
maskImg = false(size(img,1), size(img,2)); 

% class label (1: positive class, -1: negative class)
bbox      = [];   
nPos      = 0;           % # of a positive classes
while 1
    disp(['Do you continue? left click (yes) or right click (no)' ...
        'on anywhere on the image']);
 
    [~,~,button] = ginput(1);
    if button == 1
        
        nPos = nPos + 1;
        rect = getrect(fig);
    
        img = insertShape(img,'Rectangle',rect,'LineWidth',10, ...
                'Color','blue');
        imshow(img);
        bbox(nPos,:) = rect;
        
        maskImg((rect(2):(rect(2)+rect(4)-1)), ...
            (rect(1):(rect(1)+rect(3)-1))) = 1;
        
        disp(['# of Positive classes = ' int2str(nPos)]);
    else
        break;
    end
end

% random seeds
nNeg    = nPos*nNegMul;
yRange  = (1+negWinSize(end)/2):(size(img,1)-negWinSize(end)/2);
xRange  = (1+negWinSize(end)/2):(size(img,2)-negWinSize(end)/2);

iNeg    = 1;
while 1
    if (iNeg <= nNeg)
        % randomly select a center point of a negative window
        xc = xRange(randperm(length(xRange),1));
        yc = yRange(randperm(length(yRange),1));
        w  = negWinSize(randperm(numel(negWinSize),1));
        h  = negWinSize(randperm(numel(negWinSize),1));

        if sum(any(maskImg(yc-h/2:yc+h/2, xc-w/2:xc+w/2)))==0
            bbox(nPos+iNeg,:) = fix([xc-w/2 yc-h/2 w h]);
            iNeg = iNeg+1;
        end
    else
        break
    end
end
img = insertShape(img,'Rectangle',bbox(nPos+1:end,:), 'LineWidth',10, ...
        'Color','red');
imshow(img);

[~,name,~] = fileparts(orgTrImgs{idx});

imwrite(img,fullfile(folderImgTrData, ...
    [name '_posnegbox_' int2str(idx) '.jpg']));

classLabel = ones(nPos+nNeg,1);
classLabel(nPos+1:end,1) = -1;

disp('Finish the projcess !!');

save(fullfile(folderImgTrData,['DataLabel' int2str(idx)]), ...
    'classLabel','bbox');