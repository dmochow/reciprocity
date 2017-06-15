function stats = computeStatsOnly(I,head,indsActive)
% STATS=COMPUTESTATSONLY(I,HEAD,INDSACTIVE)
% compute statistics on the electric field produced by tDCS
%
% stats: structure containing the values of the x-, y-, z-, and magnitude of 
%       the electric field produced by reciprocal tDCS, as well as its normal component,
%       "focal radius" (radius of
%       sphere containing half of the total electric field), intensity of
%       targeted vertices, the normal component of the intensity of the targeted
%       vertices, and the coordinates of the peak electric field
% 
% I: applied tDCS currents (nElectrodes-by-1)
%
% head: structure that must contain lead field matrix "R" (nElectrodes x
%     (nVerticesx3)), coordinates of vertices "locs" (nVertices-by-3), surface
%     normals "surfNormals" (nVertices-by-3), mesh vertices "Vertices", mesh
%     faces "Faces", the path to an EEGLAB topoplot compatible electrode
%     location file "locFilename", and (optional) a cortical surface that can be generated
%     by the software BrainStorm "cortex"
%
% indsActive: indices of the vertices that are neurally activated

stats.E=head.R'*I;

stats.Ex=stats.E(1:3:end);
stats.Ey=stats.E(2:3:end);
stats.Ez=stats.E(3:3:end);

stats.Emag=sqrt(stats.Ex.^2+stats.Ey.^2+stats.Ez.^2);
stats.Edir=stats.Ex.*head.surfNormals(:,1)+stats.Ey.*head.surfNormals(:,2)+stats.Ez.*head.surfNormals(:,3);

%% compute focality and intensity statistics
centroid=mean( head.locs(indsActive,:) );
dists=sqrt(sum(( head.locs-repmat(centroid,size(head.locs,1),1) ).^2,2));
[nodeDistances,sortInd]=sort(dists,'ascend');

totalField=sum(stats.Emag);

directivity= cumsum( stats.Emag(sortInd) ) / totalField;
tmp=find(directivity>0.5);
if ~isempty(tmp)
    focalRadius=nodeDistances(tmp(1));
else
    focalRadius=inf;
end

intensity=mean(stats.Emag(indsActive));
intensityDir=mean(stats.Edir(indsActive));


% determine location of peak field
[~,maxInd]=max(stats.Emag);
locPeak=head.locs(maxInd,:);

stats.focalRadius=focalRadius;
stats.intensity=intensity;
stats.intensityDir=intensityDir;
stats.locPeak=locPeak;
