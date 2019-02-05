

function [x,observationVec,measurementMat,AA,BB] = experimetGeneraton(numNode,numVariable,numObservation,noiseVariance)
x = randn(numVariable,1);
observationVec = zeros(numObservation,numNode);
measurementMat = randn(numObservation,numVariable,numNode);
AA = zeros(numObservation*numNode,numVariable);
BB = zeros(numObservation*numNode,1);
for ii = 1:numNode    
    observationVec(:,ii) = measurementMat(:,:,ii)*x + sqrt(noiseVariance)*randn(numObservation,1);
    AA((ii-1)*3 + 1 : 3*ii,:) = measurementMat(:,:,ii);
    BB((ii-1)*3 + 1 : 3*ii,:) = observationVec(:,ii);
end
