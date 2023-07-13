%% Set paths and directories
% update on Jul 13 2023 by C Smid for Z Li

% Step 1: Navigate to the corticalthickness folder on your Desktop from within Matlab 
% add surfstat toolbox
addpath('surfstat')

% other functions
addpath('export_fig')

% directories for data
fsDir = 'fsaverage5_folder/'; % fsaverage5
% T0 data
T0dataDir = 'FS_registered_Clean_T0/'; % the T0 data

%% Load surfaces for fsaverage5

% Load pial and WM surfaces
lh_pial = SurfStatReadSurf([fsDir '/surf/lh.pial']);
rh_pial = SurfStatReadSurf([fsDir '/surf/rh.pial']);
pial = SurfStatAvSurf({[fsDir '/surf/lh.pial'] [fsDir '/surf/rh.pial']});
wm = SurfStatAvSurf({[fsDir '/surf/lh.white'] [fsDir '/surf/rh.white']});
inflate = SurfStatAvSurf({[fsDir '/surf/lh.inflated']  [fsDir '/surf/rh.inflated']});



% Surface: make midsurface
surf.coord    = (pial.coord+wm.coord)./2;
surf.tri      = pial.tri;
nFS = length(surf.coord); % number of vertices on the surface

SurfStatView(surf)

mask2 = SurfStatMaskCut( surf );

%% Load subject data and mask

% decision making scores and demographics
subjList = readtable('MBMF_CorticalThicknessAnalysis_Data.csv');

% general
nSubj = length(subjList.ID);
IDs = subjList.ID;%str2double(subjList.ID);
age = subjList.Age_Frac_Imp;
age_scaled = subjList.Age_Frac_Imp_z;
sex = subjList.Gender;

% 6 parameters (1 w)
it = subjList.it_P6;
lr = subjList.lr_P6;
eg = subjList.eg_P6;
w = subjList.w_lo;
st = subjList.w_hi;
repst = str2double(subjList.repst);

% seven parameters (2 ws)
%it7 = subjList.it_P6;
%lr7 = subjList.lr_P6;
%eg7 = subjList.eg_P6;
w_lo = subjList.w_lo;
w_hi = subjList.w_hi;
%st7 = subjList.st_P6;
%repst7 = subjList.repst_P6;
w_diff = subjList.w_diff;

% groups and splits
w_split = subjList.Meta_split;
w_split= num2str(w_split);
w_split(w_split=='1')='H';
w_split(w_split=='0')='L';

% EF measures
stroop = subjList.Stroop_t;
flanker = subjList.FlankerInhib_t;

flanker = flanker * -1;

% Loop over subject IDs and load cortical thickness data into large matrix
CT = zeros(nSubj, nFS);
for ii = 1:nSubj
    
    id_str = num2str(subjList.ID(ii));
    %id_str = char(subjList.ID(ii));
    %id_str = cellfun(@num2str,subjList.ID(ii),'UniformOutput',false);
    
    if length(id_str) == 1
        sub_st = 'sub-00';
    elseif length(id_str) == 2
        sub_st = 'sub-0';
    else
        sub_st = 'sub-';
    end
    
    disp(sub_st)
    disp(id_str)
    
    lh = SurfStatReadData([dataDir, sub_st, id_str, '_space-fsaverage5_desc-lh_thickness_10mm.mgh']);
    rh = SurfStatReadData([dataDir, sub_st, id_str, '_space-fsaverage5_desc-rh_thickness_10mm.mgh']);
    CT(ii,:) = [lh,rh];
end

% UPDATE 16/08/2021: from speaking to Jessica, the mask from the Surfstat
% tutorial looked the best
load '/data/mask.mat'


%% other things to report

meanthicksubj = mean( double (CT(:,mask) ),2);

% for age
f= figure;
SurfStatPlot( age, meanthicksubj, [], [], [],[], 'MarkerSize', 6 );
% set title
ylabel('Average cortical thickness')
xlabel('Age (years)')
title('Age and mean cortical thickness')
xlim([6,13])

% edit figure
set(gca,'FontSize',14)
set(gcf,'color','w')
set(gca,'Units','normalized')

%titleHandle = get(gca,'Title');
%pos = get(titleHandle,'position');
%pos1 = pos + [0 0.01 0];
%set(titleHandle,'position',pos1);

%pos = get(gca,'position');
%pos(2)=0.9*pos(2);
%set(gca,'position',pos);

h = get(gca, 'Children');
set(h(1),'Marker','o','MarkerEdgeColor','black');
set(h(2),'Color','red','linewidth',2);
set(gca,'Children',flipud(get(gca,'Children')));


export_fig  'MBSample_Thick_decrease_ages_9jun2022.jpeg' -native -m2 

% for gender
f=figure;
SurfStatPlot( sex, meanthicksubj, [], [], [], [], 'MarkerSize', 6  );
% set title
ylabel('Average cortical thickness')
xlabel('')
title('Sex and mean cortical thickness');

% edit figure
set(gca,'FontSize',14)
set(gcf,'color','w')
set(gca,'Units','normalized')

%titleHandle = get(gca,'Title');
%pos = get(titleHandle,'position');
%pos1 = pos + [0 0.01 0];
%set(titleHandle,'position',pos1);

%pos = get(gca,'position');
%pos(2)=0.9*pos(2);
%set(gca,'position',pos);

h = get(gca, 'Children');
set(h(1),'Marker','o','MarkerEdgeColor','black');
set(h(2),'Color','red','linewidth',2);
set(gca,'Children',flipud(get(gca,'Children')));

export_fig  'MBSample_Thick_Sex_9jun2022.jpeg' -native -m2 

% for w_diff - metacontrol is significantly correlated to greater mean
% thickness
f=figure,SurfStatPlot( w_diff, meanthicksubj, [], [], 'LineWidth', 2, 'MarkerSize', 12  );title('Metacontrol and mean cortical thickness');

export_fig  'MBSample_Thick_increase_metacontrol.jpeg' -native -m2 

% for w_diff - metacontrol is significantly correlated to greater mean
% thickness
f=figure,SurfStatPlot( w_diff, meanthicksubj, age, [], 'LineWidth', 2, 'MarkerSize', 12  );title('Metacontrol and mean cortical thickness (age-adjusted)');

export_fig  'MBSample_Thick_increase_metacontrol_adj_age.jpeg' -native -m2 

% for w overall
f=figure,SurfStatPlot( w, meanthicksubj, [], [], 'LineWidth', 2, 'MarkerSize', 12  );title('MB overall and mean cortical thickness');

export_fig  'MBSample_Thick_NS_MB_overall.jpeg' -native -m2

%% Linear model: age only

ageterm = term(age);

% build a model
M = 1 + ageterm;

% estimate the model parameters
slm = SurfStatLinMod(CT, M, surf);


% positive metacontrol on cortical thickness (accounting for sex and age)
slm = SurfStatT(slm, -age);

% multiple comparison correction: random field theory
[pval, peak, clus, clusid] = SurfStatP(slm, mask,0.01); % at minimum should be 0.025, now used 0.01. If put at 0.001, nothing survives
figure, SurfStatView(pval, surf, 'RFT-FWE (cap at 0.01)'); 

export_fig  'Age_effect_wholebrain_22Jun2022.jpeg' -native -m2 

% view coordinates for the cluster and peaks
term(clus) + term( SurfStatInd2Coord( clus.clusid, surf)',{'x','y','z'})
term(peak) + term( SurfStatInd2Coord( peak.vertid, surf)',{'x','y','z'})


%% Linear model: W-DIFF MODEL {Metacontrol} . wdiff + age + sex
% This is the whole-brain corrected effect of metacontrol (wdiff) on
% cortical thickness. This is our main finding from this analysis (two
% clusters survive correction)

wdiffterm = term(w_diff);
sexterm = term(sex);
ageterm = term(age);
ageterm_z = term(age_scaled);

% build a model
M = 1 + sexterm + ageterm_z + wdiffterm;

% estimate the model parameters
slm = SurfStatLinMod(CT, M, surf);


% positive metacontrol on cortical thickness (accounting for sex and age)
slm = SurfStatT(slm, w_diff);

% tvalues
tval = slm.t;
tval(~mask) = 0;
figure, SurfStatViewData(tval, surf, 'T (43 df) for metacontrol (age+sex controlled)'); SurfStatColLim([-3 3]); %colormap(hsv);

export_fig  'Wdiff_clusters_tvals_8Jun2022.jpeg' -native -m2 

% no correction
p = 1-tcdf(slm.t,slm.df);
f = figure; SurfStatViewData(p, surf, 'p-value (uncorrected)'); SurfStatColLim([0 0.05])

% multiple comparison correction: random field theory
[pval, peak, clus, clusid] = SurfStatP(slm, mask, 0.01); % at minimum should be 0.025, now used 0.01. If put at 0.001, nothing survives
figure, SurfStatView(pval, surf, 'RFT-FWE (cap at 0.01)'); 

export_fig  'Wdiff_clusters_pvals_corr_8Jun2022.jpeg' -native -m2 

% view coordinates for the cluster and peaks
term(clus) + term( SurfStatInd2Coord( clus.clusid, surf)',{'x','y','z'})
term(peak) + term( SurfStatInd2Coord( peak.vertid, surf)',{'x','y','z'})

% UPDATE 17/11/2021: i also ran this with the curvature data just to see, but there were no
% effects there.

%% check what happens in the clusters here

Meta_group = term(w_split);

% first cluster (upper)
cluster1_mask = clusid == clusid( 2452 );

% build a model
M = 1 + sexterm + ageterm_z + wdiffterm;
slm = SurfStatLinMod(CT, M, surf);
reselsROI = SurfStatResels(slm,cluster1_mask);
stat_threshold( reselsROI, sum(cluster1_mask), 1, slm.df)

YROI = mean( CT(:, cluster1_mask), 2 );

SurfStatPlot( age, YROI, Meta_group );

%% Linear model: W-DIFF MODEL {Metacontrol} . wdiff + age + sex PLUS STROOP

stroopterm = term(stroop);
wdiffterm = term(w_diff);
sexterm = term(sex);
ageterm = term(age);
ageterm_z = term(age_scaled);

% build a model
M = 1 + sexterm + ageterm_z + stroopterm + wdiffterm;

% estimate the model parameters
slm = SurfStatLinMod(CT, M, surf);

% positive metacontrol on cortical thickness (accounting for sex and age)
slm = SurfStatT(slm, w_diff);


% tvalues
tval = slm.t;
tval(~mask) = 0;
figure, SurfStatViewData(tval, surf, 'T (43 df) for metacontrol (age+sex controlled)'); SurfStatColLim([-3 3]); %colormap(hsv);

% no correction
p = 1-tcdf(slm.t,slm.df);
f = figure; SurfStatViewData(p, surf, 'p-value (uncorrected)'); SurfStatColLim([0 0.05])

% multiple comparison correction: random field theory
[pval, peak, clus, clusid] = SurfStatP(slm, mask,0.01); % at minimum should be 0.025, now used 0.01. If put at 0.001, nothing survives
figure, SurfStatView(pval, surf, 'RFT-FWE (cap at 0.01)'); 

export_fig  'Wdiff_Stroop_clusters_pvals_corr_13Jun2022.jpeg' -native -m2 

% view coordinates for the cluster and peaks
term(clus) + term( SurfStatInd2Coord( clus.clusid, surf)',{'x','y','z'})
term(peak) + term( SurfStatInd2Coord( peak.vertid, surf)',{'x','y','z'})

%% Linear model: W-DIFF MODEL {Metacontrol} . wdiff + age + sex PLUS FLANKER
% This is the whole-brain corrected effect of metacontrol (wdiff) on
% cortical thickness. This is our main finding from this analysis (two
% clusters survive correction)

flankerterm = term(flanker);
wdiffterm = term(w_diff);
sexterm = term(sex);
ageterm = term(age);
ageterm_z = term(age_scaled);

% build a model
M = 1 + sexterm + ageterm_z + flankerterm + wdiffterm;

% estimate the model parameters
slm = SurfStatLinMod(CT, M, surf);

% positive metacontrol on cortical thickness (accounting for sex and age)
slm = SurfStatT(slm, w_diff);

% tvalues
tval = slm.t;
tval(~mask) = 0;
figure, SurfStatViewData(tval, surf, 'T (43 df) for metacontrol (age+sex controlled)'); SurfStatColLim([-3 3]); %colormap(hsv);

% no correction
p = 1-tcdf(slm.t,slm.df);
f = figure; SurfStatViewData(p, surf, 'p-value (uncorrected)'); SurfStatColLim([0 0.05])

% multiple comparison correction: random field theory
[pval, peak, clus, clusid] = SurfStatP(slm, mask,0.01); % at minimum should be 0.025, now used 0.01. If put at 0.001, nothing survives
figure, SurfStatView(pval, surf, 'RFT-FWE (cap at 0.01)'); 

export_fig  'Wdiff_Flanker_clusters_pvals_corr_13Jun2022.jpeg' -native -m2 

% view coordinates for the cluster and peaks
term(clus) + term( SurfStatInd2Coord( clus.clusid, surf)',{'x','y','z'})
term(peak) + term( SurfStatInd2Coord( peak.vertid, surf)',{'x','y','z'})

%% Linear model: W-DIFF MODEL {Metacontrol} . wdiff + age + sex PLUS INHIBITION
% This is the whole-brain corrected effect of metacontrol (wdiff) on
% cortical thickness. This is our main finding from this analysis (two
% clusters survive correction)

flankerterm = term(flanker);
stroopterm = term(stroop);
wdiffterm = term(w_diff);
sexterm = term(sex);
ageterm = term(age);
ageterm_z = term(age_scaled);

% build a model
M = 1 + sexterm + ageterm_z + stroopterm + flankerterm + wdiffterm;

% estimate the model parameters
slm = SurfStatLinMod(CT, M, surf);

% positive metacontrol on cortical thickness (accounting for sex and age)
slm = SurfStatT(slm, w_diff);


% tvalues
tval = slm.t;
tval(~mask) = 0;
figure, SurfStatViewData(tval, surf, 'T (43 df) for metacontrol (age+sex controlled)'); SurfStatColLim([-3 3]); %colormap(hsv);

% no correction
p = 1-tcdf(slm.t,slm.df);
f = figure; SurfStatViewData(p, surf, 'p-value (uncorrected)'); SurfStatColLim([0 0.05])

% multiple comparison correction: random field theory
[pval, peak, clus, clusid] = SurfStatP(slm, mask,0.01); % at minimum should be 0.025, now used 0.01. If put at 0.001, nothing survives
figure, SurfStatView(pval, surf, 'RFT-FWE (cap at 0.01)'); 

% set title
title('Model-based decision making and left DLPFC')

export_fig  'Wdiff_Flanker_Stroop_clusters_pvals_corr_13Jun2022.jpeg' -native -m2 

% view coordinates for the cluster and peaks
term(clus) + term( SurfStatInd2Coord( clus.clusid, surf)',{'x','y','z'})
term(peak) + term( SurfStatInd2Coord( peak.vertid, surf)',{'x','y','z'})

%% MODEL COMPARISON for null models and models of interest

wdiffterm = term(w_diff);
sexterm = term(sex);
ageterm = term(age);
flankerterm = term(flanker);
stroopterm = term(stroop);
ageterm_z = term(age_scaled);

% Null model to WDIFF - SIGNIFICANT
% Wdiff model slightly better
slm0 = SurfStatLinMod(CT, 1 + sexterm + ageterm_z, surf);
slm1 = SurfStatLinMod(CT, 1 + sexterm + ageterm_z + wdiffterm, surf);
slm = SurfStatF(slm1,slm0);
[h,p] = ttest(slm0.SSE, slm1.SSE)
MSE0 = mean([slm0.SSE])/slm0.df
MSE1 = mean([slm1.SSE])/slm1.df
RMSE0 = (MSE0)/0.5
RMSE1 = (MSE1)/0.5
SurfStatView(slm.t.*mask, surf)

%%% COMPARE NULL AND STROOP MODEL
% No significant difference
slm0 = SurfStatLinMod(CT, 1 + sexterm + ageterm_z, surf);
slm1 = SurfStatLinMod(CT, 1 + sexterm + ageterm_z + stroopterm, surf);
slm = SurfStatF(slm1,slm0);
[h,p] = ttest(slm0.SSE, slm1.SSE)
MSE0 = mean([slm0.SSE])/slm0.df
MSE1 = mean([slm1.SSE])/slm1.df
RMSE0 = (MSE0)/0.5
RMSE1 = (MSE1)/0.5
SurfStatView(slm.t.*mask, surf)


%%% COMPARE NULL AND FLANKER MODEL
% Significant, flanker model better - actually not really any differences
slm0 = SurfStatLinMod(CT, 1 + sexterm + ageterm_z, surf);
slm1 = SurfStatLinMod(CT, 1 + sexterm + ageterm_z + flankerterm, surf);
slm = SurfStatF(slm1,slm0);
[h,p] = ttest(slm0.SSE, slm1.SSE)
MSE0 = mean([slm0.SSE])/slm0.df
MSE1 = mean([slm1.SSE])/slm1.df
RMSE0 = (MSE0)/0.5
RMSE1 = (MSE1)/0.5
SurfStatView(slm.t.*mask, surf)


%% Model comparisons of wdiff and other models with inhibition 

%%% COMPARE WDIFF AND FLANKER MODEL
% null model fits better
slm0 = SurfStatLinMod(CT, 1 + sexterm + ageterm_z + wdiffterm, surf);
slm1 = SurfStatLinMod(CT, 1 + sexterm + ageterm_z + flankerterm + wdiffterm, surf);
slm = SurfStatF(slm1,slm0);
[h,p,ci,stats] = ttest(slm0.SSE, slm1.SSE)
MSE0 = mean([slm0.SSE])/slm0.df
MSE1 = mean([slm1.SSE])/slm1.df
RMSE0 = (MSE0)/0.5
RMSE1 = (MSE1)/0.5
SurfStatView(slm.t.*mask, surf)

[pval, peak, clus, clusid] = SurfStatP(slm, mask,0.025);
term(clus) + term( SurfStatInd2Coord( clus.clusid, surf)',{'x','y','z'})
term(peak) + term( SurfStatInd2Coord( peak.vertid, surf)',{'x','y','z'})

%%% COMPARE WDIFF AND STROOP MODEL
% The full model fit slightly better
slm0 = SurfStatLinMod(CT, 1 + sexterm + ageterm_z + wdiffterm, surf);
slm1 = SurfStatLinMod(CT, 1 + sexterm + ageterm_z + stroopterm + wdiffterm, surf);
slm = SurfStatF(slm1,slm0);
[h,p,ci,stats] = ttest(slm0.SSE, slm1.SSE)
MSE0 = mean([slm0.SSE])/slm0.df
MSE1 = mean([slm1.SSE])/slm1.df
RMSE0 = (MSE0)/0.5
RMSE1 = (MSE1)/0.5
SurfStatView(slm.t.*mask, surf)

[pval, peak, clus, clusid] = SurfStatP(slm, mask,0.025);
term(clus) + term( SurfStatInd2Coord( clus.clusid, surf)',{'x','y','z'})
term(peak) + term( SurfStatInd2Coord( peak.vertid, surf)',{'x','y','z'})
 
%%% COMPARE WDIFF AND FLANKER + STROOP MODEL
% Null model # no significant differences? null model fit better
slm0 = SurfStatLinMod(CT, 1 + sexterm + ageterm_z + wdiffterm, surf);
slm1 = SurfStatLinMod(CT, 1 + sexterm + ageterm_z + flankerterm + stroopterm + wdiffterm, surf);
slm = SurfStatF(slm1,slm0);
[h,p,ci,stats] = ttest(slm0.SSE, slm1.SSE)
MSE0 = mean([slm0.SSE])/slm0.df
MSE1 = mean([slm1.SSE])/slm1.df
RMSE0 = (MSE0)/0.5
RMSE1 = (MSE1)/0.5
SurfStatView(slm.t.*mask, surf)

[pval, peak, clus, clusid] = SurfStatP(slm, mask,0.025);
term(clus) + term( SurfStatInd2Coord( clus.clusid, surf)',{'x','y','z'})
term(peak) + term( SurfStatInd2Coord( peak.vertid, surf)',{'x','y','z'})





%% Load atlases

%% Load parcellation from DK atlas
% I am using freesurfer's desikan-killiany atlas which subdivides the
% cortex into a total of 72 anatomically-defined ROIs (36 left and right).

parc = 'aparc';

% Read annotation for specified parcellaion on fsa5
[vert_lh, label_lh, ctb_lh] = read_annotation([fsDir, '/label/lh.', parc,'.annot']);
[vert_rh, label_rh, ctb_rh] = read_annotation([fsDir, '/label/rh.', parc,'.annot']);

parcFS = zeros(1,length(label_lh));
for ii = 1:length(label_lh)
    parcFS(ii) = find(ctb_lh.table(:,5) == label_lh(ii));
end
for ii = 1:length(label_rh)
    parcFS(ii+length(label_lh)) = ...
        find(ctb_rh.table(:,5) == label_rh(ii)) + length(ctb_rh.table(:,5));
end

%% load parcellation from Destrieux Atlas

parc = 'aparc.a2009s';

% Read annotation for specified parcellaion on fsa5
[vert_lh_dx, label_lh_dx, ctb_lh_dx] = read_annotation([fsDir, '/label/lh.', parc,'.annot']);
[vert_rh_dx, label_rh_dx, ctb_rh_dx] = read_annotation([fsDir, '/label/rh.', parc,'.annot']);

parcFS_dx = zeros(1,length(label_lh_dx));
for ii = 1:length(label_lh_dx)
    parcFS_dx(ii) = find(ctb_lh_dx.table(:,5) == label_lh_dx(ii));
end
for ii = 1:length(label_rh_dx)
    parcFS_dx(ii+length(label_lh_dx)) = ...
        find(ctb_rh_dx.table(:,5) == label_rh_dx(ii)) + length(ctb_rh_dx.table(:,5));
end

%% Load YEO parcellation - 7 regions
% the Yeo atlas has 16 areas in total (8 for each hemisphere, as one is
% the FreeSurfer defined medialwall)

% Read annotation for specified parcellation on fsa5
[vert_lh_ys, label_lh_ys, ctb_lh_ys] = read_annotation([fsDir, '/label/lh.Yeo2011_7Networks_N1000','.annot']);
[vert_rh_ys, label_rh_ys, ctb_rh_ys] = read_annotation([fsDir, '/label/rh.Yeo2011_7Networks_N1000','.annot']);

yeosFS = zeros(1,length(label_lh_ys));
for ii = 1:length(label_lh_ys)
    % for every vertice, we are going to find the label?
    yeosFS(ii) = find(ctb_lh_ys.table(:,5) == label_lh_ys(ii));
end
for ii = 1:length(label_rh_ys)
    % same here, appending the right numbers to this
    yeosFS(ii+length(label_lh_ys)) = find(ctb_rh_ys.table(:,5) == label_rh_ys(ii)) + length(ctb_rh_ys.table(:,5));
end


%% Load YEO parcellation - 17 regions
% the Yeo atlas has 36 areas in total (18 for each hemisphere, as one is
% the FreeSurfer defined medialwall)

% Read annotation for specified parcellation on fsa5
[vert_lh_y, label_lh_y, ctb_lh_y] = read_annotation([fsDir, '/label/lh.Yeo2011_17Networks_N1000','.annot']);
[vert_rh_y, label_rh_y, ctb_rh_y] = read_annotation([fsDir, '/label/rh.Yeo2011_17Networks_N1000','.annot']);

yeoFS = zeros(1,length(label_lh_y));
for ii = 1:length(label_lh_y)
    % for every vertice, we are going to find the label?
    yeoFS(ii) = find(ctb_lh_y.table(:,5) == label_lh_y(ii));
end
for ii = 1:length(label_rh_y)
    % same here, appending the right numbers to this
    yeoFS(ii+length(label_lh_y)) = find(ctb_rh_y.table(:,5) == label_rh_y(ii)) + length(ctb_rh_y.table(:,5));
end

%% Link cluster to brain areas

% To figure out what brain areas might be included in the clusters, I made
% an ROI from a cluster. 

% I determined the cluster ID by using the Data Cursor tool, and clicking on
% the significant clusters in the plot. (Following the ROI analysis heading
% on the SurfStat page). 

% When I do this for the cluster near the parahippocampal gyrus, I get
% several values, e.g. id: 9912, id: 605, id: 2452.

% When I sum this vector, it is always 90 - so this should be the same
% cluster. 

%% align mask ROI and parcFS - cluster 1

% areas included are , 90 voxels across four areas were included.
% The significant voxels were located in the fusiform gyrus (56.67%), 
% the entorhinal cortex (25.56%), the inferior temporal gyrus (13.33%), 
% and in the parahippocampal gyrus (4.44%). 

% ok new finding: 80 voxels across

% first cluster (lower)
maskROI = clusid == clusid( 2452 );
sum(maskROI) % this gives 90

ROIfindings(1,:) = parcFS;
ROIfindings(2,:) = maskROI;

% find all relevant areas
a = find(ROIfindings(2,:) == 1); % this matches, it's 80, which is the same as what it says in term(clus) (n of verts)

% this gives me a list of 90 areas 
locs1 = [];
for i = 1:length(a)
    % find vertex with a significant tvalue
    v_id = a(i);
    locs1(i).idx = v_id;
    locs1(i).parcFSarea = ROIfindings(1,v_id);
    locs1(i).clusid = ROIfindings(2,v_id);
    %check if left or right hemisphere (parcFS appended right to left)
    
    match = ROIfindings(1,v_id);
    if match < 37
        locs1(i).hemi = 'left';
        locs1(i).name = ctb_lh.struct_names(parcFS(v_id));
    else
        rh_id = match - 36;
        locs1(i).hemi = 'right';
        locs1(i).name = ctb_lh.struct_names(rh_id);
    end
end



%% align mask ROI and parcFS - cluster 2

% for the second cluster, located in the superior parietal cortex, the
% cluster ideas obtained via the data cursor were, e.g.: 16780, 17770,
% 14891

% 75 voxels across two areas were included. The majority of voxels belonged 
% to the right postcentral gyrus, 78.67%, and 21.33% to the superior parietal cortex. 

% first cluster (upper)
maskROI2 = clusid == clusid( 13234 );
sum(maskROI2) % this gives 75

ROIfindings2(1,:) = parcFS;
ROIfindings2(2,:) = maskROI2;

% find all relevant areas
a = find(ROIfindings2(2,:) == 1); % this matches, it's 75, which is the same as what it says in term(clus) (n of verts)

locs2 = [];
for i = 1:length(a)
    % find vertex with a significant tvalue
    v_id = a(i);
    locs2(i).idx = v_id;
    locs2(i).parcFSarea = ROIfindings2(1,v_id);
    locs2(i).clusid = ROIfindings2(2,v_id);
    %check if left or right hemisphere (parcFS appended right to left)
    
    match = ROIfindings2(1,v_id);
    if match < 37
        locs2(i).hemi = 'left';
        locs2(i).name = ctb_lh.struct_names(parcFS(v_id));
    else
        rh_id = match - 36;
        locs2(i).hemi = 'right';
        locs2(i).name = ctb_lh.struct_names(rh_id);
    end
end




%% align mask ROI and yeoFS 17 networks - cluster 1

% The paper by Shinn et al. 2015, has all the regions mapped
% https://www.frontiersin.org/files/Articles/121129/fnhum-09-00134-HTML/image_m/fnhum-09-00134-g001.jpg

% this cluster contains networks 5 (sum is 4), 9 (sum is 82) and 15 (sum is
% 4)
% network 5 is (dark green) - Dorsal attention A (posterior temporal
% occipital, superior parietal, inferior parietal occipital)

% network 9 is (light green/cream) Limbic networks (temporal pole,
% orbitofrontal)

% network 15 is... (dark blue) (default C, parahippocampal complex, ventral
% inferior parietal)

% first cluster (lower)
maskROI1_y = clusid == clusid( 2452 );
sum(maskROI1_y) % also 90 voxels

ROIfindings_y(1,:) = yeoFS;
ROIfindings_y(2,:) = maskROI1_y;

% find all relevant areas
a = find(ROIfindings_y(2,:) == 1); % this matches, it's 90

locs1_y = [];
for i = 1:length(a)
    % find vertex with a significant tvalue
    v_id = a(i);
    locs1_y(i).idx = v_id;
    locs1_y(i).yeoFSarea = ROIfindings_y(1,v_id);
    locs1_y(i).clusid = ROIfindings_y(2,v_id);
    %check if left or right hemisphere (parcFS appended right to left)
    
    match = ROIfindings_y(1,v_id);
    if match < 19
        locs1_y(i).hemi = 'left';
        locs1_y(i).name = ctb_lh_y.struct_names(yeoFS(v_id));
    else
        rh_id = match - 18;
        locs1_y(i).hemi = 'right';
        locs1_y(i).name = ctb_lh_y.struct_names(rh_id);
    end
end

% this cluster now contains network 9 and 15, so no more network 5

sum([locs1_y.yeoFSarea] == 10) % 78 is limbic network
sum([locs1_y.yeoFSarea] == 16) % 2


%% align mask ROI and yeoFS 17 networks - cluster 2


% this cluster contains networks 3 (sum is 65) and 6 (sum is 10)
% network 3 is (blue/powder blue) (somato-motor A, central sulcus secondary
% somatosensory)

% network 6 is  (dark green) (dorsal attention B, posterior temporal,
% postcentral gyrus, frontal eye fields, precentral ventral frontal)

% first cluster (upper)
maskROI2_y = clusid == clusid( 13234 );
sum(maskROI2_y)

ROIfindings_y(1,:) = yeoFS;
ROIfindings_y(2,:) = maskROI2_y;

% find all relevant areas
a = find(ROIfindings_y(2,:) == 1); % this matches, it's 130, which is the same as what it says in term(clus) (n of verts)

locs2_y = [];
for i = 1:length(a)
    % find vertex with a significant tvalue
    v_id = a(i);
    locs2_y(i).idx = v_id;
    locs2_y(i).yeoFSarea = ROIfindings_y(1,v_id);
    locs2_y(i).clusid = ROIfindings_y(2,v_id);
    %check if left or right hemisphere (parcFS appended right to left)
    
    match = ROIfindings_y(1,v_id);
    if match < 19
        locs2_y(i).hemi = 'left';
        locs2_y(i).name = ctb_lh_y.struct_names(yeoFS(v_id));
    else
        rh_id = match - 18;
        locs2_y(i).hemi = 'right';
        locs2_y(i).name = ctb_lh_y.struct_names(rh_id);
    end
end


sum([locs2_y.yeoFSarea] == 22) % 95
sum([locs2_y.yeoFSarea] == 25) % 10


%% align mask ROI and parcFS_dx - cluster 1

% first cluster (upper)
maskROI1_dx = clusid == clusid( 2452 );
sum(maskROI1_dx) % this gives 75

ROIfindings2(1,:) = parcFS_dx;
ROIfindings2(2,:) = maskROI1_dx;

% find all relevant areas
a = find(ROIfindings2(2,:) == 1); % this matches, it's 75, which is the same as what it says in term(clus) (n of verts)

locs1_dx = [];
for i = 1:length(a)
    % find vertex with a significant tvalue
    v_id = a(i);
    locs1_dx(i).idx = v_id;
    locs1_dx(i).parcFSarea = ROIfindings2(1,v_id);
    locs1_dx(i).clusid = ROIfindings2(2,v_id);
    %check if left or right hemisphere (parcFS appended right to left)
    
    match = ROIfindings2(1,v_id);
    if match < 77
        locs1_dx(i).hemi = 'left';
        locs1_dx(i).name = ctb_lh_dx.struct_names(parcFS_dx(v_id));
    else
        rh_id = match - 76;
        locs1_dx(i).hemi = 'right';
        locs1_dx(i).name = ctb_lh_dx.struct_names(rh_id);
    end
end

%% align mask ROI and parcFS_dx - cluster 2

% for the second cluster, located in the superior parietal cortex, the
% cluster ideas obtained via the data cursor were, e.g.: 16780, 17770,
% 14891

% first cluster (upper)
maskROI2_dx = clusid == clusid( 16780 );
sum(maskROI2_dx) % this gives 75

ROIfindings2(1,:) = parcFS_dx;
ROIfindings2(2,:) = maskROI2_dx;

% find all relevant areas
a = find(ROIfindings2(2,:) == 1); % this matches, it's 75, which is the same as what it says in term(clus) (n of verts)

locs2_dx = [];
for i = 1:length(a)
    % find vertex with a significant tvalue
    v_id = a(i);
    locs2_dx(i).idx = v_id;
    locs2_dx(i).parcFSarea = ROIfindings2(1,v_id);
    locs2_dx(i).clusid = ROIfindings2(2,v_id);
    %check if left or right hemisphere (parcFS appended right to left)
    
    match = ROIfindings2(1,v_id);
    if match < 77
        locs2_dx(i).hemi = 'left';
        locs2_dx(i).name = ctb_lh_dx.struct_names(parcFS_dx(v_id));
    else
        rh_id = match - 76;
        locs2_dx(i).hemi = 'right';
        locs2_dx(i).name = ctb_lh_dx.struct_names(rh_id);
    end
end


%% align mask ROI and yeoFS 7 networks - cluster 1

% so, I think that these are the networks:
% frontoparietal  = Network 6, (orange color)

% for cluster 1 for 7 networks, this is included:
% network 1 (sum is 4), network 5 (count is 86). 
% network 1 is the visual network (purple)
% network 5 is the limbic system? I think (cream/light green)

% first cluster (upper)
maskROI1_ys = clusid == clusid( 2452 );
sum(maskROI1_ys) % also 90 voxels

ROIfindings_y(1,:) = yeosFS;
ROIfindings_y(2,:) = maskROI1_ys;

% find all relevant areas
a = find(ROIfindings_y(2,:) == 1); % this matches, it's 90

locs1_ys = [];
for i = 1:length(a)
    % find vertex with a significant tvalue
    v_id = a(i);
    locs1_ys(i).idx = v_id;
    locs1_ys(i).yeoFSarea = ROIfindings_y(1,v_id);
    locs1_ys(i).clusid = ROIfindings_y(2,v_id);
    %check if left or right hemisphere (parcFS appended right to left)
    
    match = ROIfindings_y(1,v_id);
    if match < 9
        locs1_ys(i).hemi = 'left';
        locs1_ys(i).name = ctb_lh_ys.struct_names(yeosFS(v_id));
    else
        rh_id = match - 8;
        locs1_ys(i).hemi = 'right';
        locs1_ys(i).name = ctb_lh_ys.struct_names(rh_id);
    end
end

%% align mask ROI and yeoFS 7 networks - cluster 2

% for cluster 1 for 7 networks, this is included:
% network 2 (11 cause it's right) (sum is 74), network 3 (12 cause it's right) (count is 1). 
% network 2 is the somatomotor network (blue)
% network 3 is the dorsal attention network I think (dark green)

% first cluster (upper)
maskROI2_ys = clusid == clusid( 16780 );
sum(maskROI2_ys)

ROIfindings_y(1,:) = yeosFS;
ROIfindings_y(2,:) = maskROI2_ys;

% find all relevant areas
a = find(ROIfindings_y(2,:) == 1); % this matches, it's 130, which is the same as what it says in term(clus) (n of verts)

locs2_ys = [];
for i = 1:length(a)
    % find vertex with a significant tvalue
    v_id = a(i);
    locs2_ys(i).idx = v_id;
    locs2_ys(i).yeoFSarea = ROIfindings_y(1,v_id);
    locs2_ys(i).clusid = ROIfindings_y(2,v_id);
    %check if left or right hemisphere (parcFS appended right to left)
    
    match = ROIfindings_y(1,v_id);
    if match < 9
        locs2_ys(i).hemi = 'left';
        locs2_ys(i).name = ctb_lh_ys.struct_names(yeosFS(v_id));
    else
        rh_id = match - 8;
        locs2_ys(i).hemi = 'right';
        locs2_ys(i).name = ctb_lh_ys.struct_names(rh_id);
    end
end

%% ROI analysis -- testing specific brain areas

% Our approach here is to see if specific areas correlate with model-based
% decision-making and/or metacontrol. 

% With Boris's help, I am now using the residual cortical thickness
% controlled for age in the models. 

%% extracting means of ROIs

% From meeting with Boris 13/09/21:
% build a linear model, but do not fit it to the mean thickness vector of that
% area. no multiple comparisons. fontol parietal network - define based on
% Yeo : mean thickness of each cluster, and at the end a model for each of
% the clusters.

% linear model is to use the mean score, w and age and sex, 
% means will be a number of regions per converted p-values

% yeo networks in yeo paper of 2011 

% QUESTION: Should I not stick to the same atlas for all these analyses?

%%% 28/09/21 BORIS: you can plot an area on the surface
ROI_both = double(parcFS == 15 | parcFS == 51);
SurfStatViewData(ROI_both, surf)

% left only
SurfStatViewData(double(parcFS == 15), surf)

%%% 28/09/21 Boris: this approach shows residual cortical thickness controlled for
% age, correlate this with metacontrol
slm1 = SurfStatLinMod(CT, 1 + ageterm);

% residual cortical thickness
CTC = CT - slm1.X*slm1.coef;

% THIS IS RESIDUAL THICKNESS CONTROLLED FOR AGE
meanCTC = mean(CTC(:, idxOFL),2);

% this actually controls for age --> This is your ROI analysis, controlled
% for age
figure,SurfStatPlot(w_diff, meanCTC, [], [], 'LineWidth', 2, 'MarkerSize', 12)


%% FRONTAL AREAS

%% dorso-lateral prefrontal cortex (rostral middle frontal cortex in DK)

% these are areas 28 (left) and 64 (right)

%%% 28/09/21 BORIS: you can plot an area on the surface? --> to check if i get the region right
ROI_both = double(parcFS == 28 | parcFS == 64);

% manually set colors
map = [0.5 0.5 0.5
    0.0297 0.7082 0.8163];

map = [0.827 0.827 0.827
    0.5044 0.7993 0.3480];

SurfStatViewData(ROI_both, surf, 'DLPFC'); colormap(map);

export_fig  'DLPFC_doubleROI.jpeg' -native -m2 


% left only
SurfStatViewData(double(parcFS == 28), surf)

%%% left DLPFC
idxOFL = find(parcFS == 28);

meanCT = mean(CT(:, idxOFL),2);

%%% 28/09/21 Boris: this approach shows residual cortical thickness controlled for
% age, correlate this with metacontrol
slm1 = SurfStatLinMod(CT, 1 + ageterm_z);

% residual cortical thickness
CTC = CT - slm1.X*slm1.coef;

% THIS IS RESIDUAL THICKNESS CONTROLLED FOR AGE
meanCTC = mean(CTC(:, idxOFL),2);

% overall model-based DM
figure,SurfStatPlot(w, meanCTC, [], [], 'MarkerSize', 6)

% set title
ylabel('Residual Cortical Thickness')
xlabel('T = -1.52, df = 42, p = 0.135')
title('Model-based decision making and left DLPFC')

% edit figure
set(gca,'FontSize',14)
set(gcf,'color','w')
set(gca,'Units','normalized')

%titleHandle = get(gca,'Title');
%pos = get(titleHandle,'position');
%pos1 = pos + [0 0.01 0];
%set(titleHandle,'position',pos1);

%pos = get(gca,'position');
%pos(2)=0.9*pos(2);
%set(gca,'position',pos);

h = get(gca, 'Children');
set(h(1),'Marker','o','MarkerEdgeColor','black');
set(h(2),'Color','red','linewidth',2);
set(gca,'Children',flipud(get(gca,'Children')));

export_fig  'Scatter_lDLPFC_w_8jun2022.jpeg' -native -m3 

% metacontrol
figure,SurfStatPlot(w_diff, meanCTC, [], [], 'MarkerSize', 6)

% set title
ylabel('Residual Cortical Thickness')
xlabel('T = 2.61, df = 42, p = 0.012')
title('Metacontrol and left DLPFC')

% edit figure
set(gca,'FontSize',14)
set(gcf,'color','w')
set(gca,'Units','normalized')

%titleHandle = get(gca,'Title');
%pos = get(titleHandle,'position');
%pos1 = pos + [0 0.01 0];
%set(titleHandle,'position',pos1);

%pos = get(gca,'position');
%pos(2)=0.9*pos(2);
%set(gca,'position',pos);

h = get(gca, 'Children');
set(h(1),'Marker','o','MarkerEdgeColor','black');
set(h(2),'Color','red','linewidth',2);
set(gca,'Children',flipud(get(gca,'Children')));

export_fig  'Scatter_lDLPFC_wdiff_8jun2022.jpeg' -native -m3 


%%% right
SurfStatViewData(double(parcFS == 64), surf)

% right DLPFC
idxOFR = find(parcFS == 64);

meanCT = mean(CT(:, idxOFR),2);

%%% 28/09/21 Boris: this approach shows residual cortical thickness controlled for
% age, correlate this with metacontrol
slm1 = SurfStatLinMod(CT, 1 + ageterm_z);

% residual cortical thickness
CTC = CT - slm1.X*slm1.coef;

% THIS IS RESIDUAL THICKNESS CONTROLLED FOR AGE
meanCTC = mean(CTC(:, idxOFR),2);

% model-based DM
figure,SurfStatPlot(w, meanCTC, [], [], 'MarkerSize', 6)

% set title
ylabel('Residual Cortical Thickness')
xlabel('T = -1.71, df = 42, p = 0.094')
title('Model-based decision making and right DLPFC')

% edit figure
set(gca,'FontSize',14)
set(gcf,'color','w')
set(gca,'Units','normalized')

%titleHandle = get(gca,'Title');
%pos = get(titleHandle,'position');
%pos1 = pos + [0 0.01 0];
%set(titleHandle,'position',pos1);

%pos = get(gca,'position');
%pos(2)=0.9*pos(2);
%set(gca,'position',pos);

h = get(gca, 'Children');
set(h(1),'Marker','o','MarkerEdgeColor','black');
set(h(2),'Color','red','linewidth',2);
set(gca,'Children',flipud(get(gca,'Children')));

export_fig  'Scatter_rDLPFC_w_8jun2022.jpeg' -native -m3 


% metacontrol
figure,SurfStatPlot(w_diff, meanCTC, [], [], 'MarkerSize', 6)

% set title
ylabel('Residual Cortical Thickness')
xlabel('T = 3.00, df = 42, p = 0.005')
title('Metacontrol and right DLPFC')

% edit figure
set(gca,'FontSize',14)
set(gcf,'color','w')
set(gca,'Units','normalized')

%titleHandle = get(gca,'Title');
%pos = get(titleHandle,'position');
%pos1 = pos + [0 0.01 0];
%set(titleHandle,'position',pos1);

%pos = get(gca,'position');
%pos(2)=0.9*pos(2);
%set(gca,'position',pos);

h = get(gca, 'Children');
set(h(1),'Marker','o','MarkerEdgeColor','black');
set(h(2),'Color','red','linewidth',2);
set(gca,'Children',flipud(get(gca,'Children')));

export_fig  'Scatter_rDLPFC_wdiff_8jun2022.jpeg' -native -m3 


%% dorso-lateral prefrontal cortex (rostral middle frontal cortex in DK) for INHIBITION MEASURES

%%% NS FOR FLANKER
% these are areas 28 (left) and 64 (right)
%%% left DLPFC
idxOFL = find(parcFS == 28);

%%% 28/09/21 Boris: this approach shows residual cortical thickness controlled for
% age, correlate this with metacontrol
slm1 = SurfStatLinMod(CT, 1 + ageterm_z);

% residual cortical thickness
CTC = CT - slm1.X*slm1.coef;

% THIS IS RESIDUAL THICKNESS CONTROLLED FOR AGE
meanCTC = mean(CTC(:, idxOFL),2);

%%% FLANKER - NS
figure,SurfStatPlot(flanker, meanCTC, [], [], 'MarkerSize', 6)

% set title
ylabel('Residual Cortical Thickness')
xlabel('T = -1.74, df = 42, p = 0.089')
title('Flanker Inhibition and left DLPFC')

% edit figure
set(gca,'FontSize',14)
set(gcf,'color','w')
set(gca,'Units','normalized')

h = get(gca, 'Children');
set(h(1),'Marker','o','MarkerEdgeColor','black');
set(h(2),'Color','red','linewidth',2);
set(gca,'Children',flipud(get(gca,'Children')));

export_fig  'Scatter_lDLPFC_flanker_13jun2022.jpeg' -native -m3 

%%% STROOP - NS
figure,SurfStatPlot(stroop, meanCTC, [], [], 'MarkerSize', 6)

% set title
ylabel('Residual Cortical Thickness')
xlabel('T = -0.25, df = 42, p = 0.807')
title('Stroop Inhibition and left DLPFC')

% edit figure
set(gca,'FontSize',14)
set(gcf,'color','w')
set(gca,'Units','normalized')


h = get(gca, 'Children');
set(h(1),'Marker','o','MarkerEdgeColor','black');
set(h(2),'Color','red','linewidth',2);
set(gca,'Children',flipud(get(gca,'Children')));

export_fig  'Scatter_lDLPFC_stroop_13jun2022.jpeg' -native -m3 


% right DLPFC
idxOFR = find(parcFS == 64);

%%% 28/09/21 Boris: this approach shows residual cortical thickness controlled for
% age, correlate this with metacontrol
slm1 = SurfStatLinMod(CT, 1 + ageterm_z);

% residual cortical thickness
CTC = CT - slm1.X*slm1.coef;

% THIS IS RESIDUAL THICKNESS CONTROLLED FOR AGE
meanCTC = mean(CTC(:, idxOFR),2);

%%% FLANKER - NS BUT AT TREND
figure,SurfStatPlot(flanker, meanCTC, [], [], 'MarkerSize', 6)

% set title
ylabel('Residual Cortical Thickness')
xlabel('T = -3.01, df = 42, p = 0.004')
title('Flanker Inhibition and right DLPFC')

% edit figure
set(gca,'FontSize',14)
set(gcf,'color','w')
set(gca,'Units','normalized')

h = get(gca, 'Children');
set(h(1),'Marker','o','MarkerEdgeColor','black');
set(h(2),'Color','red','linewidth',2);
set(gca,'Children',flipud(get(gca,'Children')));

export_fig  'Scatter_rDLPFC_flanker_13jun2022.jpeg' -native -m3 


%%% STROOP - NS
figure,SurfStatPlot(stroop, meanCTC, [], [], 'MarkerSize', 6)

% set title
ylabel('Residual Cortical Thickness')
xlabel('T = 0.25, df = 42, p = 0.801')
title('Stroop Inhibition and right DLPFC')

% edit figure
set(gca,'FontSize',14)
set(gcf,'color','w')
set(gca,'Units','normalized')

h = get(gca, 'Children');
set(h(1),'Marker','o','MarkerEdgeColor','black');
set(h(2),'Color','red','linewidth',2);
set(gca,'Children',flipud(get(gca,'Children')));

export_fig  'Scatter_rDLPFC_stroop_13jun2022.jpeg' -native -m3 



%% control relationship between the DLPFC and metacontrol and flanker?

%%% NS FOR FLANKER
% these are areas 28 (left) and 64 (right)
%%% left DLPFC
idxOFL = find(parcFS == 28);

%%% 28/09/21 Boris: this approach shows residual cortical thickness controlled for
% age, correlate this with metacontrol
slm1 = SurfStatLinMod(CT, 1 + ageterm_z + flankerterm);

% residual cortical thickness
CTC = CT - slm1.X*slm1.coef;

% THIS IS RESIDUAL THICKNESS CONTROLLED FOR AGE
meanCTC = mean(CTC(:, idxOFL),2);

%%% FLANKER - NS
figure,SurfStatPlot(wdiffterm, meanCTC, [], [], 'MarkerSize', 6)

% set title
ylabel('Residual Cortical Thickness')
xlabel('T = 1.99, df = 42, p = 0.053')
title('Metacontrol controlled by Flanker and left DLPFC')

% edit figure
set(gca,'FontSize',14)
set(gcf,'color','w')
set(gca,'Units','normalized')

h = get(gca, 'Children');
set(h(1),'Marker','o','MarkerEdgeColor','black');
set(h(2),'Color','red','linewidth',2);
set(gca,'Children',flipud(get(gca,'Children')));

export_fig  'Scatter_lDLPFC_meta_flanker_13jun2022.jpeg' -native -m3 

% right DLPFC
idxOFR = find(parcFS == 64);

%%% 28/09/21 Boris: this approach shows residual cortical thickness controlled for
% age, correlate this with metacontrol
slm1 = SurfStatLinMod(CT, 1 + ageterm_z + flankerterm);

% residual cortical thickness
CTC = CT - slm1.X*slm1.coef;

% THIS IS RESIDUAL THICKNESS CONTROLLED FOR AGE
meanCTC = mean(CTC(:, idxOFR),2);

%%% FLANKER - NS BUT AT TREND
figure,SurfStatPlot(wdiffterm, meanCTC, [], [], 'MarkerSize', 6)

% set title
ylabel('Residual Cortical Thickness')
xlabel('T = 2.05, df = 42, p = 0.046')
title('Metacontrol controlled by Flanker and right DLPFC')

% edit figure
set(gca,'FontSize',14)
set(gcf,'color','w')
set(gca,'Units','normalized')

h = get(gca, 'Children');
set(h(1),'Marker','o','MarkerEdgeColor','black');
set(h(2),'Color','red','linewidth',2);
set(gca,'Children',flipud(get(gca,'Children')));

export_fig  'Scatter_rDLPFC_meta_flanker_13jun2022.jpeg' -native -m3 




%% --- STRUCTURAL COVARIANCE ANALYSIS

%% Correlating thickness of a seed region to the rest of the cortex

% UPDATE 25/11/2021: This is following Boris' paper of empathy in women



%% Linear model: W-DIFF MODEL {Metacontrol} . wdiff + age + sex
% This is the whole-brain corrected effect of metacontrol (wdiff) on
% cortical thickness. This is our main finding from this analysis (two
% clusters survive correction)

wdiffterm = term(w_diff);
sexterm = term(sex);
ageterm = term(age);

% build a model
M = 1 + sexterm + ageterm + wdiffterm;

% estimate the model parameters
slm = SurfStatLinMod(CT, M, surf);

% positive metacontrol on cortical thickness (accounting for sex and age)
slm = SurfStatT(slm, w_diff);

%find(slm.t == max(slm.t)) % 8668
%Yseed = double (CT(:,81887));

% tvalues
tval = slm.t;
tval(~mask) = 0;
figure, SurfStatViewData(tval, surf, 'T (43 df) for metacontrol (age+sex controlled)'); SurfStatColLim([-3 3]); %colormap(hsv);

% no correction
p = 1-tcdf(slm.t,slm.df);
f = figure; SurfStatViewData(p, surf, 'p-value (uncorrected)'); SurfStatColLim([0 0.05])

% multiple comparison correction: random field theory
[pval, peak, clus, clusid] = SurfStatP(slm, mask,0.01); % at minimum should be 0.025, now used 0.01. If put at 0.001, nothing survives
figure, SurfStatView(pval, surf, 'RFT-FWE (cap at 0.01)'); 

% view coordinates for the cluster and peaks
term(clus) + term( SurfStatInd2Coord( clus.clusid, surf)',{'x','y','z'})
term(peak) + term( SurfStatInd2Coord( peak.vertid, surf)',{'x','y','z'})

% UPDATE 13/12/2021: to check location of the seed I just looked at the
% value in parcFS for the id number (5485)
% id 11455 is in the corpus callosum... so this means that the children
% with more metacontrol have higher connections to the corpus callosum??

% the peak for cluster 1 is at 5485 - which is inside the entorhinal cortex
Yseed1 = double (CT(:,5485));

% the maximum peak for metacontrol in cluster 1
figure, SurfStatPlot( w_diff, Yseed1 );

% the peak for cluster 2 is at 13242 - which is in the superior parietal
% cortex
Yseed2 = double (CT(:,13242));

% the maximum peak for metacontrol in cluster 1
figure, SurfStatPlot( w_diff, Yseed2 );

%% Connectivity
% Let's take the seed as the peak of the metacontrol effect above. Is the cortical
% thickness at this point correlated with the rest of the cortical surface,
% allowing for an effect of metacontrol and age?

slm = SurfStatLinMod( CT, 1 + wdiffterm + term(Yseed1), surf );
slm = SurfStatT( slm, Yseed1 );

SurfStatView( slm.t, surf, 'Connectivity' );
SurfStatColLim( [0 2.5] );
%Of course the T statistic is infinite at the seed itself, so we needed to reset the colour limits: 

% A large part of the brain is correlated with the seed. However I would 
% claim that pure connectivity as above is of little scientific interest.
% A more interesting question is how the connectivity changes with say 
% high and low metacontrol users. To answer this, add an interaction between seed and the rest of 
% the model, then look at the interaction between seed and metacontrol groups:

group = term(wsplit_term);
%group = term(wsplit_term3);

% the peak difference in the entorhinal cortex between low and high
% metacontrol users
slm = SurfStatLinMod( CT, ( 1 + group + ageterm )*( 1 + term(Yseed1) ), surf );
slm = SurfStatT( slm, Yseed1.*group.H - Yseed1.*group.L  );
SurfStatView( slm.t.*mask, surf, 'High metacontrol - Low metacontrol connectivity' );

export_fig  'Structural_Connectivity_Seed1_groupdifference_brain.jpeg' -native -m3 

% next we find the peak of this difference.
% id 11455 is in the corpus callosum... so this means that the children
% with more metacontrol have higher connections from the entorhinal cortex to the corpus callosum??
find(slm.t == max(slm.t)) % 11455 % 18784
Yseed1_seed = double (CT(:,11455));

% LOCATION OF THIS SEED = CORPUS CALLOSUM

% MODEL

M = (1 + ageterm) * (1 + term(Yseed1));
SurfStatPlot( Yseed1, Yseed1_seed, M, group);

% there seems to be increased connectivity for this point for high and low
% metacontrol users.
% We see that there is more connectivity in high metacontrol users than low
% metacontrol users. Which is why the T statistic was positive. the F
% statistic is as expected too

figure,SurfStatPlot(Yseed1, Yseed1_seed, M, group, 'MarkerSize', 6)

% set title
ylabel('Left seed region adjusted for age')
xlabel('Seed * Metacontrol(H/L): F(1,39) = 13.48, p = 0.001')
title('High metacontrol - Low metacontrol connectivity')

% edit figure
set(gca,'FontSize',14)
set(gcf,'color','w')
set(gca,'Units','normalized')

titleHandle = get(gca,'Title');
pos = get(titleHandle,'position');
pos1 = pos + [0 0.01 0];
set(titleHandle,'position',pos1);

h = get(gca, 'Children');
set(h(1),'Marker','o');
set(h(2),'linewidth',2);
set(h(3),'Marker','o');
set(h(4),'linewidth',2);

export_fig  'Structural_Connectivity_Seed1_groupdifference.jpeg' -native -m3 


%%% SECOND SEED BASED ON PARIETAL CORTEX
slm = SurfStatLinMod( CT, 1 + wdiffterm + term(Yseed2), surf );
slm = SurfStatT( slm, Yseed2 );

SurfStatView( slm.t, surf, 'Connectivity' );
SurfStatColLim( [0 2.5] );
%Of course the T statistic is infinite at the seed itself, so we needed to reset the colour limits: 

% A large part of the brain is correlated with the seed. However I would 
% claim that pure connectivity as above is of little scientific interest.
% A more interesting question is how the connectivity changes with say 
% gender. To answer this, add an interaction between seed and the rest of 
% the model, then look at the interaction between seed and metacontrol groups:

group = term(wsplit_term);
group = term(wsplit_term3);

slm = SurfStatLinMod( CT, ( 1 + group + ageterm )*( 1 + term(Yseed2) ), surf );
slm = SurfStatT( slm, Yseed2.*group.H - Yseed2.*group.L  );
SurfStatView( slm.t.*mask, surf, 'High metacontrol - Low metacontrol connectivity' );

export_fig  'Structural_Connectivity_Seed2_groupdifference_brain.jpeg' -native -m3 

find(slm.t == max(slm.t)) % 11547 % 18781
Yseed2_seed = double (CT(:,18781));

M = (1 + ageterm) * (1 + term(Yseed2));
%SurfStatPlot( Yseed2, Yseed2_seed, M, group);

% there seems to be increased connectivity for this point for high and low
% metacontrol users

% Again we see that high metacontrol users show higher structural
% connectivity for this point compared to low metacontrol users. 

figure,SurfStatPlot(Yseed2, Yseed2_seed, M, group, 'MarkerSize', 6)

% set title
ylabel('Right seed region adjusted for age')
xlabel('Seed * Metacontrol(H/L): F(1,39) = 2.21, p = 0.145')
title('High metacontrol - Low metacontrol connectivity')

% edit figure
set(gca,'FontSize',14)
set(gcf,'color','w')
set(gca,'Units','normalized')

titleHandle = get(gca,'Title');
pos = get(titleHandle,'position');
pos1 = pos + [0 0.01 0];
set(titleHandle,'position',pos1);

h = get(gca, 'Children');
set(h(1),'Marker','o');
set(h(2),'linewidth',2);
set(h(3),'Marker','o');
set(h(4),'linewidth',2);

export_fig  'Structural_Connectivity_Seed2_groupdifference.jpeg' -native -m3 




%% REPEATING THIS BUT FOR THE MODEL INCLUDING INHIBITION

%% Correlating thickness of a seed region to the rest of the cortex

% UPDATE 25/11/2021: This is following Boris' paper of empathy in women

%% Linear model: W-DIFF MODEL {Metacontrol} . wdiff + age + sex
% This is the whole-brain corrected effect of metacontrol (wdiff) on
% cortical thickness. This is our main finding from this analysis (two
% clusters survive correction)

flankerterm = term(flanker);
stroopterm = term(stroop);
wdiffterm = term(w_diff);
sexterm = term(sex);
ageterm = term(age);

% build a model
M = 1 + sexterm + ageterm + wdiffterm + flankerterm + stroopterm;

% estimate the model parameters
slm = SurfStatLinMod(CT, M, surf);

% positive metacontrol on cortical thickness (accounting for sex and age)
slm = SurfStatT(slm, w_diff);

%find(slm.t == max(slm.t)) % 8668
%Yseed = double (CT(:,81887));

% tvalues
tval = slm.t;
tval(~mask) = 0;
figure, SurfStatViewData(tval, surf, 'T (43 df) for metacontrol (age+sex controlled)'); SurfStatColLim([-3 3]); %colormap(hsv);

% no correction
p = 1-tcdf(slm.t,slm.df);
f = figure; SurfStatViewData(p, surf, 'p-value (uncorrected)'); SurfStatColLim([0 0.05])

% multiple comparison correction: random field theory
[pval, peak, clus, clusid] = SurfStatP(slm, mask,0.01); % at minimum should be 0.025, now used 0.01. If put at 0.001, nothing survives
figure, SurfStatView(pval, surf, 'RFT-FWE (cap at 0.01)'); 

% view coordinates for the cluster and peaks
term(clus) + term( SurfStatInd2Coord( clus.clusid, surf)',{'x','y','z'})
term(peak) + term( SurfStatInd2Coord( peak.vertid, surf)',{'x','y','z'})

% the peak for cluster 1 is at 5485
Yseed1 = double (CT(:,5485));

% the maximum peak for metacontrol in cluster 1
SurfStatPlot( w_diff, Yseed1 );

% the peak for cluster 2 is at 13242
Yseed2 = double (CT(:,13242));

% the maximum peak for metacontrol in cluster 1
SurfStatPlot( w_diff, Yseed2 );

%% Connectivity
%Let's take the seed as the peak of the metacontrol effect above. Is the cortical
% thickness at this point correlated with the rest of the cortical surface,
% allowing for an effect of metacontrol and age?

slm = SurfStatLinMod( CT, 1 + wdiffterm + term(Yseed1), surf );
slm = SurfStatT( slm, Yseed1 );

SurfStatView( slm.t, surf, 'Connectivity' );
SurfStatColLim( [0 2.5] );
%Of course the T statistic is infinite at the seed itself, so we needed to reset the colour limits: 

% A large part of the brain is correlated with the seed. However I would 
% claim that pure connectivity as above is of little scientific interest.
% A more interesting question is how the connectivity changes with say 
% gender. To answer this, add an interaction between seed and the rest of 
% the model, then look at the interaction between seed and metacontrol groups:

group = term(wsplit_term);

slm = SurfStatLinMod( CT, ( 1 + group + ageterm )*( 1 + term(Yseed1) ), surf );
slm = SurfStatT( slm, Yseed1.*group.H - Yseed1.*group.L  );
SurfStatView( slm.t.*mask, surf, 'High metacontrol - Low metacontrol connectivity' );

export_fig  'Structural_Connectivity_Seed1_groupdifference_brain.jpeg' -native -m3 

find(slm.t == max(slm.t)) % 11455
Yseed1_seed = double (CT(:,11455));

M = (1 + ageterm) * (1 + term(Yseed1));
SurfStatPlot( Yseed1, Yseed1_seed, M, group);

% there seems to be increased connectivity for this point for high and low
% metacontrol users.
% We see that there is more connectivity in high metacontrol users than low
% metacontrol users. Which is why the T statistic was positive. the F
% statistic is as expected too

figure,SurfStatPlot(Yseed1, Yseed1_seed, M, group, 'MarkerSize', 6)

% set title
ylabel('Left seed region adjusted for age')
xlabel('Seed * Metacontrol(H/L): F(1,39) = 13.48, p = 0.001')
title('High versus low metacontrol users structural connectivity (left cluster)')

% edit figure
set(gca,'FontSize',14)
set(gcf,'color','w')
set(gca,'Units','normalized')

titleHandle = get(gca,'Title');
pos = get(titleHandle,'position');
pos1 = pos + [0 0.01 0];
set(titleHandle,'position',pos1);

h = get(gca, 'Children');
set(h(1),'Marker','o');
set(h(2),'linewidth',2);
set(h(3),'Marker','o');
set(h(4),'linewidth',2);

export_fig  'Structural_Connectivity_Seed1_groupdifference.jpeg' -native -m3 


%%% SECOND SEED BASED ON PARIETAL CORTEX
slm = SurfStatLinMod( CT, 1 + wdiffterm + term(Yseed2), surf );
slm = SurfStatT( slm, Yseed2 );

SurfStatView( slm.t, surf, 'Connectivity' );
SurfStatColLim( [0 2.5] );
%Of course the T statistic is infinite at the seed itself, so we needed to reset the colour limits: 

% A large part of the brain is correlated with the seed. However I would 
% claim that pure connectivity as above is of little scientific interest.
% A more interesting question is how the connectivity changes with say 
% gender. To answer this, add an interaction between seed and the rest of 
% the model, then look at the interaction between seed and metacontrol groups:

group = term(wsplit_term);

slm = SurfStatLinMod( CT, ( 1 + group + ageterm )*( 1 + term(Yseed2) ), surf );
slm = SurfStatT( slm, Yseed2.*group.H - Yseed2.*group.L  );
SurfStatView( slm.t.*mask, surf, 'High metacontrol - Low metacontrol connectivity' );

export_fig  'Structural_Connectivity_Seed2_groupdifference_brain.jpeg' -native -m3 

find(slm.t == max(slm.t)) % 11547
Yseed2_seed = double (CT(:,11455));

M = (1 + ageterm) * (1 + term(Yseed2));
%SurfStatPlot( Yseed2, Yseed2_seed, M, group);

% there seems to be increased connectivity for this point for high and low
% metacontrol users

% Again we see that high metacontrol users show higher structural
% connectivity for this point compared to low metacontrol users. 

figure,SurfStatPlot(Yseed2, Yseed2_seed, M, group, 'MarkerSize', 6)

% set title
ylabel('Right seed region adjusted for age')
xlabel('Seed * Metacontrol(H/L): F(1,39) = 2.21, p = 0.145')
title('High versus low metacontrol users structural connectivity (right cluster)')

% edit figure
set(gca,'FontSize',14)
set(gcf,'color','w')
set(gca,'Units','normalized')

titleHandle = get(gca,'Title');
pos = get(titleHandle,'position');
pos1 = pos + [0 0.01 0];
set(titleHandle,'position',pos1);

h = get(gca, 'Children');
set(h(1),'Marker','o');
set(h(2),'linewidth',2);
set(h(3),'Marker','o');
set(h(4),'linewidth',2);

export_fig  'Structural_Connectivity_Seed2_groupdifference.jpeg' -native -m3 


%% DLPFC MEAN CORTICAL THICKNESS AND STRUCTURAL CONNECTIVITY

% right DLPFC
idxOFR = find(parcFS == 64);

%%% 28/09/21 Boris: this approach shows residual cortical thickness controlled for
% age, correlate this with metacontrol
slm1 = SurfStatLinMod(CT, 1 + ageterm);

% residual cortical thickness
CTC = CT - slm1.X*slm1.coef;

% THIS IS RESIDUAL THICKNESS CONTROLLED FOR AGE
meanCTC = mean(CTC(:, idxOFR),2);

%%% 
figure,SurfStatPlot(w_diff, meanCTC, [], [], 'MarkerSize', 6)

% include this in a model?
DLPFC_term = term(meanCTC);

group = term(wsplit_term);

slm = SurfStatLinMod( CT, ( 1 + group )*( 1 + term(DLPFC_term) ), surf );
slm = SurfStatT( slm, meanCTC.*group.H - meanCTC.*group.L  );
SurfStatView( slm.t.*mask, surf, 'High metacontrol - Low metacontrol connectivity' );

% tvalues
tval = slm.t;
tval(~mask) = 0;
figure, SurfStatViewData(tval, surf, 'T (43 df) for metacontrol (age+sex controlled)'); SurfStatColLim([-3 3]); %colormap(hsv);

% no correction
p = 1-tcdf(slm.t,slm.df);
f = figure; SurfStatViewData(p, surf, 'p-value (uncorrected)'); SurfStatColLim([0 0.05])

% multiple comparison correction: random field theory
[pval, peak, clus, clusid] = SurfStatP(slm, mask,0.025); % at minimum should be 0.025, now used 0.01. If put at 0.001, nothing survives
figure, SurfStatView(pval, surf, 'RFT-FWE (cap at 0.025)'); 

% view coordinates for the cluster and peaks
term(clus) + term( SurfStatInd2Coord( clus.clusid, surf)',{'x','y','z'})
term(peak) + term( SurfStatInd2Coord( peak.vertid, surf)',{'x','y','z'})

find(slm.t == max(slm.t)) % 11455
Yseed_DLPFC = double (CT(:,11455));

M = (1 + ageterm) * (1 + term(Yseed1));
SurfStatPlot( meanCTC, Yseed_DLPFC, M, group);


% or just looking at the model?
M = 1 + wdiffterm + DLPFC_term;
slm = SurfStatLinMod(CT, M, surf)
slm = SurfStatT(slm, meanCTC);


%% Linear model: W-DIFF MODEL {Metacontrol} . wdiff + age + sex
% This is the whole-brain corrected effect of metacontrol (wdiff) on
% cortical thickness. This is our main finding from this analysis (two
% clusters survive correction)

lrterm = term(lr);
sexterm = term(sex);
ageterm_z = term(age_scaled);

% build a model
M = 1 + sexterm + ageterm_z + lrterm;

% estimate the model parameters
slm = SurfStatLinMod(CT, M, surf);


% positive metacontrol on cortical thickness (accounting for sex and age)
slm = SurfStatT(slm, lr);

% tvalues
tval = slm.t;
tval(~mask) = 0;
figure, SurfStatViewData(tval, surf, 'T (42 df) for lr (age+sex controlled)'); SurfStatColLim([-3 3]); %colormap(hsv);

% no correction
p = 1-tcdf(slm.t,slm.df);
f = figure; SurfStatViewData(p, surf, 'p-value (uncorrected)'); SurfStatColLim([0 0.05])

% multiple comparison correction: random field theory
[pval, peak, clus, clusid] = SurfStatP(slm, mask,0.01); % at minimum should be 0.025, now used 0.01. If put at 0.001, nothing survives
figure, SurfStatView(pval, surf, 'RFT-FWE (cap at 0.01)'); 

export_fig  'lr_clusters_pvals_corr_8Jun2022.jpeg' -native -m2 

% other parameters
itterm = term(it);

% build a model
M = 1 + sexterm + ageterm_z + itterm;

% estimate the model parameters
slm = SurfStatLinMod(CT, M, surf);

% positive metacontrol on cortical thickness (accounting for sex and age)
slm = SurfStatT(slm, it);

% tvalues
tval = slm.t;
tval(~mask) = 0;
figure, SurfStatViewData(tval, surf, 'T (42 df) for lr (age+sex controlled)'); SurfStatColLim([-3 3]); %colormap(hsv);

export_fig  'it_clusters_tvals_8Jun2022.jpeg' -native -m2 

% no correction
p = 1-tcdf(slm.t,slm.df);
f = figure; SurfStatViewData(p, surf, 'p-value (uncorrected)'); SurfStatColLim([0 0.05])

% multiple comparison correction: random field theory
[pval, peak, clus, clusid] = SurfStatP(slm, mask,0.01); % at minimum should be 0.025, now used 0.01. If put at 0.001, nothing survives
figure, SurfStatView(pval, surf, 'RFT-FWE (cap at 0.01)'); 

% other parameters
egterm = term(eg);

% build a model
M = 1 + sexterm + ageterm_z + egterm;

% estimate the model parameters
slm = SurfStatLinMod(CT, M, surf);

% positive metacontrol on cortical thickness (accounting for sex and age)
slm = SurfStatT(slm, eg);

% tvalues
tval = slm.t;
tval(~mask) = 0;
figure, SurfStatViewData(tval, surf, 'T (42 df) for lr (age+sex controlled)'); SurfStatColLim([-3 3]); %colormap(hsv);

% no correction
p = 1-tcdf(slm.t,slm.df);
f = figure; SurfStatViewData(p, surf, 'p-value (uncorrected)'); SurfStatColLim([0 0.05])

% multiple comparison correction: random field theory
[pval, peak, clus, clusid] = SurfStatP(slm, mask,0.01); % at minimum should be 0.025, now used 0.01. If put at 0.001, nothing survives
figure, SurfStatView(pval, surf, 'RFT-FWE (cap at 0.01)'); 

export_fig  'eg_clusters_pvals_8Jun2022.jpeg' -native -m2 

% other parameters
stterm = term(st);

% build a model
M = 1 + sexterm + ageterm_z + stterm;

% estimate the model parameters
slm = SurfStatLinMod(CT, M, surf);

% positive metacontrol on cortical thickness (accounting for sex and age)
slm = SurfStatT(slm, st);

% tvalues
tval = slm.t;
tval(~mask) = 0;
figure, SurfStatViewData(tval, surf, 'T (42 df) for lr (age+sex controlled)'); SurfStatColLim([-3 3]); %colormap(hsv);

% no correction
p = 1-tcdf(slm.t,slm.df);
f = figure; SurfStatViewData(p, surf, 'p-value (uncorrected)'); SurfStatColLim([0 0.05])

% multiple comparison correction: random field theory
[pval, peak, clus, clusid] = SurfStatP(slm, mask,0.01); % at minimum should be 0.025, now used 0.01. If put at 0.001, nothing survives
figure, SurfStatView(pval, surf, 'RFT-FWE (cap at 0.01)'); 

export_fig  'st_clusters_pvals_8Jun2022.jpeg' -native -m2 
