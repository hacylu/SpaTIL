function [ features,featureNames ] = getSpaTILFeatures( lympCentroids,nonLympCentroids,alpha )
%GETSPATILFEATURES Gets a set of spatial quantitative features (SpaTIL) 
% relating to the spatial architecture of TILs, co-localization of TILs and 
% cancer nuclei, and density variation of TIL clusters.

if nargin<3
    alpha=.49;
end
r=.185;

lympNodes = struct('centroid_r',lympCentroids(:,2)','centroid_c',lympCentroids(:,1)');
nonLympNodes = struct('centroid_r',nonLympCentroids(:,2)','centroid_c',nonLympCentroids(:,1)');

[~,~,~,~,lympMatrix,~,~] = construct_nodesCluster(lympNodes, alpha, r);
[~,~,~,~,nonLympMatrix,~,~] = construct_nodesCluster(nonLympNodes, alpha, r);

[~,~,lympMembers] = networkComponents(lympMatrix);
[~,~,nonLympMembers] = networkComponents(nonLympMatrix);

  [ feat5,featDescription5 ] = getClusterAreaMeasures(lympMembers,lympNodes,nonLympMembers,nonLympNodes);
  [ feat6,featDescription6 ] = getClusterDensityMeasures(lympMembers,lympNodes,nonLympMembers,nonLympNodes);
  [ feat7,featDescription7 ] = getIntersectedArea( lympMembers,lympNodes,nonLympMembers,nonLympNodes );
  
  [ feat8_1,featDescription81 ]  = extractNNearestClusterProperties(lympMembers,lympNodes,nonLympMembers,nonLympNodes,1);
  [ feat8_2,featDescription82 ]  = extractNNearestClusterProperties( lympMembers,lympNodes,nonLympMembers,nonLympNodes,2);
  [ feat8_3,featDescription83 ]  = extractNNearestClusterProperties( lympMembers,lympNodes,nonLympMembers,nonLympNodes,3);
  [ feat8_4,featDescription84 ]  = extractNNearestClusterProperties( lympMembers,lympNodes,nonLympMembers,nonLympNodes,4);
  [ feat8_5,featDescription85 ]  = extractNNearestClusterProperties( lympMembers,lympNodes,nonLympMembers,nonLympNodes,5);
  feat8 = [feat8_1,feat8_2,feat8_3,feat8_4,feat8_5];
  featDescription8 = [featDescription81,featDescription82,featDescription83,featDescription84,featDescription85];

LymPoints = {};
LymConcaveCenter = [];
numLymMembers=length(lympMembers);

for j = 1:numLymMembers
    a = lympMembers{j};
    LymConcaveC = lympNodes.centroid_c(a)';
    LymConcaveR = lympNodes.centroid_r(a)';
    LymPoints(j) = {[LymConcaveC,LymConcaveR]};
    ConcaveCenter = [mean(LymConcaveC),mean(LymConcaveR)];
    LymConcaveCenter = [LymConcaveCenter;ConcaveCenter];
end

NucleiPoints = {};
NucleiConcaveCenter = [];

numNucleiMembers=length(nonLympMembers);
for j = 1:numNucleiMembers
    a = nonLympMembers{j};
    NucleiConcaveC = nonLympNodes.centroid_c(a)';
    NucleiConcaveR = nonLympNodes.centroid_r(a)';
    NucleiPoints(j) = {[NucleiConcaveC,NucleiConcaveR]};
    ConcaveCenter = [mean(NucleiConcaveC),mean(NucleiConcaveR)];
    NucleiConcaveCenter = [NucleiConcaveCenter;ConcaveCenter];
end

[feat1,featDescription1] = extractLocalClusterDensityRelationship(LymPoints,NucleiPoints,LymConcaveCenter,NucleiConcaveCenter);
[feat2,featDescription2] = extractLocalMinimalLymClusterEncompassed(LymPoints,NucleiPoints,LymConcaveCenter,NucleiConcaveCenter);
[feat3,featDescription3] = extracRatioOfMinimalLymClusterEncompassed(LymPoints,NucleiPoints,LymConcaveCenter,NucleiConcaveCenter);

%I=imread('D:\German\Data\TMA_UH\TMA1\123_1.png');
%[NucleiedgesConnection,NuceliClusterCenterKept,LymedgesConnection,LymClusterCenterKept] = extractClusterGraph_visual(LymPoints,NucleiPoints,LymConcaveCenter,NucleiConcaveCenter,I);
[NucleiedgesConnection,NuceliClusterCenterKept,LymedgesConnection,LymClusterCenterKept] = extractClusterGraph(LymPoints,NucleiPoints,LymConcaveCenter,NucleiConcaveCenter);


% a = logical(NucleiedgesConnection);
% a = tril(a)'|a|triu(a)';
% NucleiedgesConnection = a;

NucleiG = graph(NucleiedgesConnection);
% a = logical(LymedgesConnection);
% a = tril(a)'|a|triu(a)';
% LymedgesConnection = a;
LymG = graph(LymedgesConnection);
NucleiImportance = computeIntegralNodesImportance(NucleiG,NuceliClusterCenterKept);
LymImportance = computeIntegralNodesImportance(LymG,LymClusterCenterKept);
[feat4,featDescription4] = findIndexOfNearestLym(NucleiImportance,NuceliClusterCenterKept,LymImportance,LymClusterCenterKept);

features = [feat1,feat2,feat3,feat4,feat5,feat6,feat7,feat8];
featureNames = [featDescription1,featDescription2,featDescription3,featDescription4,featDescription5,featDescription6,featDescription7,featDescription8];

end

