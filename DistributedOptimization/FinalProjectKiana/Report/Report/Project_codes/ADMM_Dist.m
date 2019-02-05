function [errorNormalizedDist] = ADMM_Dist(adjMatrixOrig,rho,numNode,numberVariable,numIteration,A,b,xHatCen)

xHatMatrixNew = zeros(numberVariable,numNode);
xHatMatrixNew2 = zeros(numberVariable,numNode);
alphaMat = zeros(numberVariable,numNode);
errorNormalizedDist = zeros(numIteration,1);
% xHatMatrixIterDist = zeros(numberVariable,numNode,numIteration);
lastRCVDMat = zeros(numberVariable,numNode,numNode);

for jj = 1:numIteration
   
        for ii = 1:numNode
            xHatMatrixOld = xHatMatrixNew;
            nodeInd = ii;
            nodeIndex = nodeInd;
            numNeighbor = sum(adjMatrixOrig(nodeInd,:));
            lastRCVDMat(:,:,nodeInd) = repmat(adjMatrixOrig(nodeInd,:),size(xHatMatrixNew,1),1).*xHatMatrixOld;
            alphaMat(:,nodeInd) = alphaMat(:,nodeInd) + rho*(numNeighbor*xHatMatrixOld(:,nodeInd) - sum(lastRCVDMat(:,:,nodeInd),2)) ;
            ACV = A(:,:,nodeIndex);
            bCV = b(:,nodeIndex);
            xHatMatrixNew2(:,nodeIndex) = inv(ACV'*ACV + 2*rho*numNeighbor*eye(numberVariable))*( ACV'*bCV - alphaMat(:,nodeIndex) + ...
            rho*(numNeighbor*xHatMatrixOld(:,nodeIndex) + sum(lastRCVDMat(:,:,nodeIndex),2)));
        end
        xHatMatrixNew = xHatMatrixNew2;
       
    errorNormalizedDist(jj,:) = sqrt(sum(sum((xHatMatrixNew - repmat(xHatCen,1,numNode)).^2)))/sqrt(sum(sum(repmat(xHatCen,1,numNode).^2)));  
    
end




