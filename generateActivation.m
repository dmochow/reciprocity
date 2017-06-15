function [V,o] = generateActivation(X,N,S,head,sigma) 
%   V=GENERATEACTIVATION(X,N,S,HEAD,SIGMA)
%   simulate activation patterns of seed points in X, sized by N, scaled by
%   S, using the head model head and noise standard deviation sigma
%
%   X: nSeeds x 3 matrix of seed positions
%   N: nSeeds x 1 vector of number of vertices in each seed
%   S: nSeeds x 1 vector of intensity (factor of 1e-9 m*A)
%   head: head model structure: must have fields "locs" (nVertices-by-3), "R" (nElectrodes-by-(nVerticesx3)) , and
%       "surfNormals" (nVertices-by-3)
%   sigma: standard deviation of white noise added to electrodes
%
%   V: nSeeds x 1: the simulated scalp potentials
%   o: output struct containing x-, y-, and z- components of cortical
%   current density, as well as its magnitude

nSeedPoints=size(X,1);

J=zeros(3*size(head.Vertices,1),1);
Jo=1e-9; % 1 nM*A

for s=1:nSeedPoints
   
    xSeed=X(s,:);
    theseDists=sqrt(sum(( head.locs-repmat(xSeed,size(head.locs,1),1) ).^2,2));
    [~,sortInd]=sort(theseDists,'ascend');
    inds=sortInd(1:N(s));
    
    J(3*inds-2)=S(s)*Jo*head.surfNormals(inds,1); % x-direction
    J(3*inds-1)=S(s)*Jo*head.surfNormals(inds,2); % y-direction
    J(3*inds)=S(s)*Jo*head.surfNormals(inds,3); % z-direction
    
end

sigma=sigma*std(head.R*J);
W=sigma*randn(size(head.R,1),1);
V=head.R*J+W;

o.J=J;
o.Jx=J(1:3:end);
o.Jy=J(2:3:end);
o.Jz=J(3:3:end);
o.Jmag=sqrt(o.Jx.^2+o.Jy.^2+o.Jz.^2);



