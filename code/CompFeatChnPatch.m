function trFeat = CompFeatChnPatch(itgPatchSet,haarPatch, nHaar)

nChannel    = size(itgPatchSet,3);
matrixSize  = [size(itgPatchSet,1) size(itgPatchSet,2)]; 
featureIdx  = 1;
trFeat      = zeros(1,nChannel*nHaar);
for ii = 1 : nChannel
    patch   = itgPatchSet(:,:,ii);
    for jj=1:nHaar
        posMat  = haarPatch{jj,1}(:, 1:2);
        operMat = haarPatch{jj,1}(:, 3);
        
        idx = sub2ind(matrixSize,posMat(:,2),posMat(:,1));
        trFeat(1,featureIdx)  = patch(idx)'*operMat;
        featureIdx  = featureIdx + 1;
    end
end
