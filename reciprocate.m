function I = reciprocate(V,head,Imax)
%I=RECIPROCATE(R,V,IMAX) determine the TES montage that targets the sources expressed in
%scalp potentials V, generated as a linear combination of the columns of R,
%and constrained to a maximum delivered current of Imax
%   R: channels x (locs*3) lead-field matrix
%   V: channels x 1 vector of scalp potentials
%   Imax: total current delivered in A (defaults to 0.002) 

if nargin<3, Imax=0.002; end

RR=head.R*(head.R)';

rankRR=rank(RR);
if rankRR<size(head.R,1),
    warning('R is rank-deficient.  Results will likely be wrong');
    % try a regularized inverse
    warning('Trying to compute a regularized inverse...');
    %nElectrodes=size(RR,1);
    %lambda=1e0;
    %I = inv( RR+lambda*eye(nElectrodes) ) \ V;
    I = regInv( RR,100 ) * V;
else
    I=RR\V;  % R is well-conditioned
end


Itotal=sum(abs([I; sum(I)]));
I=I/(Itotal/(2*Imax));

% %% compute the electric field from reciprocal TES
% o.E=head.R'*I;
% o.Ex=o.E(1:3:end);
% o.Ey=o.E(2:3:end);
% o.Ez=o.E(3:3:end);
% o.Emag=sqrt(o.Ex.^2+o.Ey.^2+o.Ez.^2);
% o.Edir=o.Ex.*head.surfNormals(:,1)+o.Ey.*head.surfNormals(:,2)+o.Ez.*head.surfNormals(:,3);
% 
% %% compute focality and intensity statistics
% indsActive=find(Jmag);
% centroid=mean( head.locs(indsActive,:) ); 
% nActive=numel(indsActive);
% if ~isempty(Jmag)
%     
%     dists=sqrt(sum(( head.locs-repmat(centroid,size(head.locs,1),1) ).^2,2));
%     [nodeDistances,sortInd]=sort(dists,'ascend');
%     
%     totalField=sum(o.Emag);
%     
%     directivity= cumsum( o.Emag(sortInd) ) / totalField;
%     tmp=find(directivity>0.5);
%     if ~isempty(tmp)
%         focalRadius=nodeDistances(tmp(1));
%     else
%         focalRadius=inf;
%     end
%     
%     intensity=mean(o.Emag(indsActive));
%     intensityDir=mean(o.Edir(indsActive));
% 
%     o.focalRadius=focalRadius;
%     o.intensity=intensity;
%     o.intensityDir=intensityDir;
% end

% opts: to be implemented in the future (L1)

end

