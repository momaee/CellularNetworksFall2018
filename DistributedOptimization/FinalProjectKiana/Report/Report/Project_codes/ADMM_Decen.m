function [ errorNormalizedDecen ] = ADMM_Decen( c,numAgent,numVariable,numIteration,measurementMat,observationVec,xHatCen)
errorNormalizedDecen = zeros(numIteration,1);
xHatMatrixNew = zeros(numVariable,numAgent);

Yold = zeros(numVariable,numAgent);
Zold = zeros(numVariable,numAgent);

for jj = 1:numIteration
        for ii = 1:numAgent
            AgentIndex = ii;
            ACV = measurementMat(:,:,AgentIndex);
            bCV = observationVec(:,AgentIndex);
            xHatMatrixNew(:,AgentIndex) = ((ACV)'*((ACV)) + c*eye(numVariable))\(c*Zold(:,AgentIndex) - Yold(:,AgentIndex) + ((ACV)')*bCV);
            
        end
        
        Znew = repmat(mean(xHatMatrixNew + (1/c)*Yold,2),1,numAgent);
        Ynew = Yold + c*(xHatMatrixNew-Znew);
        
        Yold = Ynew;
        Zold = Znew; 
                        
 % xHatIterDecenter(:,:,jj) = xHatMatrixNew;           
         
    errorNormalizedDecen(jj,:) = sqrt(sum(sum((xHatMatrixNew - repmat(xHatCen,1,numAgent)).^2)))/sqrt(sum(sum(repmat(xHatCen,1,numAgent).^2)));  
    
    
    
    
    
    
end






