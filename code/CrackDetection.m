function [crackBW, flag] = CrackDetection(patchOrg, Param)

patchOrgMedian = medfilt2(patchOrg, [5 5]);
[mask, centPt, isObj] = RemoveObjArea(patchOrgMedian);

if isObj ==0; 
    flag = 0;
    crackBW   = [];
    return
end;

% crop patch and mask so that
[patch, xRange, yRange]  = MoveCenter(patchOrg, centPt);
[mask , ~]       = MoveCenter(mask, centPt);

% find Radon range
[R,~] = radon(mask, Param.theta);
radonRange  = logical(R);

STATS   = regionprops(mask, 'EquivDiameter');

% make a donut mask
mask    = bwdist(mask);
mask    = and(mask < STATS.EquivDiameter, mask > Param.borderMargin);

options = struct('FrangiScaleRange', Param.FrangiScaleRange, ...
    'FrangiScaleRatio', Param.FrangiScaleRatio, ...
    'FrangiBetaOne', Param.FrangiBetaOne, ...
    'FrangiBetaTwo', Param.FrangiBetaTwo, ...
    'verbose',Param.verbose, ...
    'BlackWhite',Param.BlackWhite);

crackPatch=FrangiFilter2D(patch*255,options);

BW = crackPatch>Param.FrangiFilterThresh;
BW = imclearborder(BW, 8);
BW = bwareaopen(BW, Param.minEdgeArea);

BW = (BW .* mask) > 0;
stats = regionprops(BW, 'Centroid', 'Area', 'Perimeter','PixelList',  ...
    'Perimeter','MajorAxisLength', 'MinorAxisLength');

col_obj = [];
row_obj = [];
AreaTmp = 0;
for index=1:length(stats)
    ratioMajMin = stats(index).MajorAxisLength/stats(index).MinorAxisLength;
    if (Param.rgd_ratioMajMin  < ratioMajMin)
        if AreaTmp ==0; AreaTmp = stats(index).Area; end;
        if stats(index).Area < AreaTmp; continue; end;
        col_obj = stats(index).PixelList(1,1);
        row_obj = stats(index).PixelList(1,2);
    end
end

BW1    = bwselect(BW,col_obj,row_obj,8);

[R,~]  = radon(BW1,Param.theta);

[~,ind]=max(R);
[~,ind1]=max(max(R));
cc = ind(ind1);
rr = ind1;

if (radonRange(cc,rr) ~= 0)
    crackBW = false(size(patchOrg,1),size(patchOrg,2));
    crackBW(yRange,xRange) = BW1;
    flag = 1;
else
    flag = 0;
    crackBW = [];
end

end
function [mask, cent, isObj]= RemoveObjArea(patch)

isObj   = 1;
img     = patch;

if nnz(img) == 0
    mask    = zeros(size(img));
    isObj   = 0;
    cent    = [ 0 0];
    return
end

edgeCanny = edge(img, 'canny', [0.2 0.5]*median(img(:)));
% figure, imshow(edgeCanny), title('binary gradient mask');

se = strel('disk',8,0);
BWsdil      = imdilate(edgeCanny, se );
% figure, imshow(BWserode), title('dilated gradient mask');

% BWsfill = imfill(BWsdil, 'holes'); figure, imshow(BWsfill),
% title('dilated gradient mask');

BWnobord = imclearborder(BWsdil, 8);
% figure, imshow(BWnobord), title('cleared border image');

LB = ceil(0.01 * size(img,1) * size(img,2));
UB = ceil(0.5  * size(img,1) * size(img,2));
Iout = xor(bwareaopen(BWnobord,LB),  bwareaopen(BWnobord,UB));
% figure, imshow(Iout), title('cleared border image');

L = bwlabel(Iout,8);

if max(L(:)) >= 2       % in case, detected objects are more than 1.
    
    range = size(L,1)/8;
    idx = L(fix(size(L,1)/2)-range:fix(size(L,1)/2)+range, ...
        fix(size(L,2)/2)-range:fix(size(L,2)/2)+range);
    
    if  nnz(idx)==0 || (length(unique(idx))>2);
        mask    = zeros(size(img));
        isObj   = 0;
        cent    = [ 0 0];
        return
    end;
    tmp = idx(find(idx));
    Iout(:,:) = 0; Iout(find(L==tmp(1))) = 1;
end

objImg  = bwconvhull(Iout);
stats   = regionprops(objImg,'Centroid');

mask = objImg;

if isempty(stats); cent = [ 0 0];   isObj = 0;
    
else cent = stats.Centroid; end;

end
function [patch, xRange, yRange] = MoveCenter(patch, centPt)

patchSize = size(patch);

ptX    = centPt(1);
ptY    = centPt(2);

if ptX < patchSize(2)-ptX;
    xRange = 1 : ptX*2;
else
    xRange = -patchSize(2)+2*ptX+1 : patchSize(2);
end

if ptY < patchSize(1)-ptY;
    yRange = 1 : ptY*2;
else
    yRange = -patchSize(1)+2*ptY+1 : patchSize(1);
end

yRange  = fix(yRange);
xRange  = fix(xRange);

patch   = patch(yRange,xRange);

end
