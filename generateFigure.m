% please report any bugs to jdmochowski AT ccny DOT cuny DOT edu


%% generate figure comparing original reciprocity to L1-constrained reciprocity

% if you find this code useful, please cite:  
% Dmochowski, J. P., Koessler, L., Norcia, A. M., Bikson, M., & Parra, L. C. (2017). 
% Optimal use of EEG recordings to target active brain areas with transcranial electrical stimulation. 
% NeuroImage.


clear all; close all; clc
addpath(genpath('.'));

% the following two lines load in the results of targeting EEG activations
% from all 15,002 possible targets in this head model
allDataFilename='./data/allLocationsResults.mat';
load(allDataFilename);

% load colormaps
load('redgreyblue.mat');
cmrgb=cm;
load('redgrey.mat');
cmrg=cm;
cmusa=jmaColors('usadarkblue');

% some parameters for visualization
plotrad=0.6;
headrad=0.545;
zoomFactor=2;
viewAngle=[-90 30];

% parameters for performing reciprocal tDCS
Imax=0.002; % maximum current = 2 mA
noiseStd=0.01; % variance of noise added to EEG sensors (regularizes covariance)
c=1e10;  % key parameter "c" for the L^1 constrained version of reciprocity

% load in the head model here
headFile='./data/icbm';
load(headFile);

%% generate activation
% two "seed" points: one in left visual cortex, one in the right
X=[-0.06964 -0.008389 0.06919; -0.07039 -0.00001825 0.06583 ];
N=[25; 25]; % for each seed point, choose the 25 nearest vertices to include in activation
S=[1 1]; % strength of each activation in nm*A
K=63; % number of electrodes
[V,o1] = generateActivation(X,N,S,head,noiseStd);

%% perform unconstrained reciprocity
Irec = reciprocate( V*c,head,Imax); % core function
indsActive=find(o1.Jmag); % determine which nodes were (neurally) activated in the head model 
statsRec = computeStatsOnly(Irec,head,indsActive); % comput statistics on the resulting electric field from tDCS

%% perform L^1 constrained reciprocity
H=2*(head.R*head.R'); % quadratic term
f=-2*V*c; % linear term
maxIter=500; % max number of iterations in quadratic programming search: increase this if not converged after 500
t=2*Imax; % upper-bound on L^1 norm of tDCS currents (4 mA for a 2 mA max total current delivered)
Il1= tibshirani(H,f,t,maxIter);
statsL1 = computeStatsOnly(Il1,head,indsActive);

%% visualize results
figure
hs(1)=subplot(2,4,1);
fCortex = patch( 'Vertices',head.Vertices, 'Faces',head.Faces,...
    'FaceVertexCData',o1.Jmag*1e9,...  %gCortex.colors,...   %
    'FaceColor','interp', 'FaceLighting','gouraud', 'BackFaceLighting','unlit', ...
    'EdgeColor','none', 'DiffuseStrength',0.7, 'SpecularStrength',0.05,...
    'SpecularExponent',5, 'SpecularColorReflectance',0.5 , ...
    'LineStyle','none');
axis off;
axis equal;
view(viewAngle);
camlight('headlight')
colormap(gca,colormap(cmrg));
zoom(zoomFactor);
hcb(1)=colorbar('south');

hs(2)=subplot(2,4,2);
topoplot(V*1e6,head.locFilename,'plotrad',plotrad,'headrad',headrad,'electrodes','off','numcontour',0);
colormap(gca,colormap(cmusa));
hcb(2)=colorbar('south');


hs(3)=subplot(2,4,3);
topoplot_dc(Irec*1e3,head.locFilename,'plotrad',plotrad,'headrad',headrad,'electrodes','off');
colormap(gca,colormap(cmusa));
hcb(3)=colorbar('south');

hs(4)=subplot(2,4,4);
fCortex = patch( 'Vertices',head.Vertices, 'Faces',head.Faces,...
    'FaceVertexCData',statsRec.Emag,...  %gCortex.colors,...   %
    'FaceColor','interp', 'FaceLighting','gouraud', 'BackFaceLighting','unlit', ...
    'EdgeColor','none', 'DiffuseStrength',0.7, 'SpecularStrength',0.05,...
    'SpecularExponent',5, 'SpecularColorReflectance',0.5 , ...
    'LineStyle','none');
axis off;
axis equal;
view(viewAngle);
camlight('headlight')
colormap(gca,colormap('jet'));
zoom(zoomFactor);
hcb(4)=colorbar('south');

% 
hs(5)=subplot(2,4,5);
topoplot_dc(Il1*1e3,head.locFilename,'plotrad',plotrad,'headrad',headrad,'electrodes','off');
colormap(gca,colormap(cmusa));
hcb(5)=colorbar('south');   

hs(6)=subplot(2,4,6);
fCortex = patch( 'Vertices',head.Vertices, 'Faces',head.Faces,...
    'FaceVertexCData',statsL1.Emag,...  %gCortex.colors,...   %
    'FaceColor','interp', 'FaceLighting','gouraud', 'BackFaceLighting','unlit', ...
    'EdgeColor','none', 'DiffuseStrength',0.7, 'SpecularStrength',0.05,...
    'SpecularExponent',5, 'SpecularColorReflectance',0.5 , ...
    'LineStyle','none');
axis off;
axis equal;
view(viewAngle);
camlight('headlight')
colormap(gca,colormap('jet'));
zoom(zoomFactor);
hcb(6)=colorbar('south');

%% show the results of targeting to all possible locations 
% these results are not computed above
% they are preloaded from disk
hs(7)=subplot(247);
hbar(1)=bar([1 2],[mean(focalRadiusTrue) mean(focalRadiusL1)]*1e2,'k'); hold on
herr(1)=errorbar([1 2],[mean(focalRadiusTrue) mean(focalRadiusL1)]*1e2,[std(focalRadiusTrue) std(focalRadiusL1)]*1e2,'k','LineStyle','none');

set(hbar(1),'FaceColor',[0.7 0.7 0.7]);
set(hbar(1),'EdgeColor','none');
ylabel('Focality (cm)');
xlim([0.25 2.75]);
axis square;
set(gca,'XTickLabel',{'Unc.','L1'});
box off

hs(8)=subplot(248);
hbar(2)=bar([1 2],[mean(intensityTrue) mean(intensityL1)],'k'); hold on
herr(2)=errorbar([1 2],[mean(intensityTrue) mean(intensityL1)],[std(intensityTrue) std(intensityL1)],'k','LineStyle','none');


% statistical tests
[pi]=signrank(intensityTrue,intensityL1);
[pf]=signrank(focalRadiusTrue,focalRadiusL1);

%% esthetics
set(hbar(2),'FaceColor',[0.7 0.7 0.7]);
set(hbar(2),'EdgeColor','none');
ylabel('Intensity (V/m)');
axis square
xlim([0.25 2.75]);
set(gca,'XTickLabel',{'Unc.','L1'});
box off;

% equate colormap limits where appropriate
% currents
clim5=get(hs(5),'CLim');
set(hs(3),'CLim',clim5);

% fields
clim6=get(hs(6),'CLim');
set(hs(4),'Clim',[0 clim6(2)]);
set(hs(6),'Clim',[0 clim6(2)]);

% now the messy part: fix positions and add colorbars

% move entire bottom row up
delUp=0.1;
for s=5:8
    pos=get(hs(s),'Position');
    set(hs(s),'Position',[pos(1) pos(2)+delUp pos(3) pos(4)]);
end

% move entire second column left
delLeft=0.025;
for s=[2 6]
    pos=get(hs(s),'Position');
    set(hs(s),'Position',[pos(1)-delLeft pos(2) pos(3) pos(4)]);
end

% move entire third column right
delRight=0.025;
for s=[3 7]
    pos=get(hs(s),'Position');
    set(hs(s),'Position',[pos(1)+delRight pos(2) pos(3) pos(4)]);
end

% create some room for bottom right subpanel y-axis labels
delLabel=0.02;

pos=get(hs(7),'Position');
set(hs(7),'Position',[pos(1)-delLabel pos(2) pos(3) pos(4)]);

pos=get(hs(8),'Position');
set(hs(8),'Position',[pos(1)+delLabel pos(2) pos(3) pos(4)]);

%
set(get(hs(2),'Title'),'String','Activation');
set(get(hs(3),'Title'),'String','Unconstrained reciprocity');
set(get(hs(5),'Title'),'String','L1-constrained reciprocity');

htit=get(hs(2),'Title');
moveTitle(htit,-0.6,0,0);

htit=get(hs(3),'Title');
moveTitle(htit,0.6,0,0);

htit=get(hs(5),'Title');
moveTitle(htit,0.6,0,0);

hcb(1) = moveHorizontalColobar( hs(1) , hcb(1) , 'nM \cdot A', 0.5, 0.1 );
hcb(4) = moveHorizontalColobar( hs(4) , hcb(4) , 'V/m', 0.5, 0.1 );
hcb(6) = moveHorizontalColobar( hs(6) , hcb(6) , 'V/m', 0.5, 0.1 );

hcb(2) = moveHorizontalColobar( hs(2) , hcb(2) , '\muV', 0.5, 0.1 );
hcb(3) = moveHorizontalColobar( hs(3) , hcb(3) , 'mA', 0.5, 0.1 );
hcb(5) = moveHorizontalColobar( hs(5) , hcb(5) , 'mA', 0.5, 0.1 );
