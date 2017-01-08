function itgPatchSet = CompIntChnImg(img)

% a set of integral images
itgPatchSet = zeros(size(img,1)+1,size(img,2)+1, 1);

% LUV transformation---------------------------------------------------
% LUV and remove L component
patchLUV = rgbConvert(img, 'luv'); 
itgPatchSet(:,:,1) = integralImage(patchLUV(:,:,2));
itgPatchSet(:,:,2) = integralImage(patchLUV(:,:,3));

% HSV transformation---------------------------------------------------
% HSV and remove V component
patchHSV = rgbConvert(img, 'hsv'); 
itgPatchSet(:,:,3) = integralImage(patchHSV(:,:,1));
itgPatchSet(:,:,4) = integralImage(patchHSV(:,:,2));

% Gradient histogram---------------------------------------------------
[M,O]       = gradientMag(single(img));
nOrients    = 6;
H           = gradientHist(M,O,1,nOrients,0,0,0); 
for ii=1:nOrients
    itgPatchSet(:,:,4+ii) = integralImage(H(:,:,ii)); 
end; 

% Gradient img --------------------------------------------------------
itgPatchSet(:,:,4+nOrients+1) = integralImage(M);
