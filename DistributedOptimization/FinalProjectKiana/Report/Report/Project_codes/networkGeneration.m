

function [adjacencyMat] = networkGeneration(treeFalg,branchFactor,regularFlag,degree,numNode,numEdge)
if treeFalg==1
    edgeList = canonicalNets(numNode,'tree',branchFactor);  % Constructing edge lists for simple canonical graphs, ex: trees and lattices.
    adjacencyMat = edgeL2adj(edgeList);
elseif regularFlag==1
    edgeList = kregular(numNode,degree);    % Create a k-regular graph.
    adjacencyMat = edgeL2adj(edgeList);
else
    connectedFlag = 0;
    while connectedFlag==0
        adjacencyMat = randomGraph(numNode,.5,numEdge);    % Random graph construction routine.
        connectedFlag = isConnected(adjacencyMat);
    end
end


