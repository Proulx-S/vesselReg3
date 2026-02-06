clear all
close all
% clc
% figure('Menu','none','ToolBar','none');

projectName = 'vesselReg3';
scriptDir   = 'doIt_vesByVesVolPrc';
%%%%%%%%%%%%%%%%%%%%%
%% Set up environment
%%%%%%%%%%%%%%%%%%%%%

% Detect computing environment
os   = char(java.lang.System.getProperty('os.name'));
host = char(java.net.InetAddress.getLocalHost.getHostName);
user = char(java.lang.System.getProperty('user.name'));

% Setup folders
if strcmp(os,'Linux') && strcmp(host,'takoyaki') && strcmp(user,'sebp')
    envId      = 1;
    storageDrive = '/local/users/Proulx-S/';
    scratchDrive = '/scratch/users/Proulx-S/';
    projectCode    = fullfile(scratchDrive, projectName);        if ~exist(projectCode,'dir');    mkdir(projectCode);    end
    projectStorage = fullfile(storageDrive, projectName);        if ~exist(projectStorage,'dir'); mkdir(projectStorage); end
    projectScratch = fullfile(scratchDrive, projectName, 'tmp'); if ~exist(projectScratch,'dir'); mkdir(projectScratch); end
    toolDir        = '/scratch/users/Proulx-S/tools';            if ~exist(toolDir,'dir');        mkdir(toolDir);        end
else
    % envId = 2;
    % storageDrive = '/Users/sebastienproulx/bassWrap-reg';
    % scratchDrive = '/Users/sebastienproulx/bassWrap-reg';
    % projectCode    = fullfile(scratchDrive, projectName);        if ~exist(projectCode,'dir');    mkdir(projectCode);    end
    % projectStorage = fullfile(storageDrive, projectName);        if ~exist(projectStorage,'dir'); mkdir(projectStorage); end
    % projectScratch = fullfile(scratchDrive, projectName, 'tmp'); if ~exist(projectScratch,'dir'); mkdir(projectScratch); end
    % toolDir        = '/Users/sebastienproulx/tools';             if ~exist(toolDir,'dir');        mkdir(toolDir);        end
end

% Load dependencies and set paths
%%% matlab util (contains my matlab git wrapper)
tool = 'util'; toolURL = 'https://github.com/Proulx-S/util.git';
if ~exist(fullfile(toolDir, tool), 'dir'); system(['git clone ' toolURL ' ' fullfile(toolDir, tool)]); end; addpath(genpath(fullfile(toolDir,tool)))
%%% matlab others
tool = 'freesurfer'; subTool = 'matlab'; repoURL = 'https://github.com/freesurfer/freesurfer.git';
gitClone(repoURL, fullfile(toolDir, tool), subTool);
%%% bassWrap-reg tools
tool = 'bassWrap-reg';
if exist(fullfile(toolDir, tool), 'dir'); addpath(fullfile(toolDir, tool)); end

%%% neurodesk
switch envId
    case 1
        global src    
        setenv('SINGULARITY_BINDPATH',strjoin({projectCode projectStorage projectScratch toolDir},','));
        %%%% vesselboost
        src.vesselboost = 'ml vesselboost/1.0.0';
        system([src.vesselboost '; prediction.py --help > /dev/null'],'-echo');
        vesselBoostModel = fullfile(projectScratch,'manual_0429');
        system([src.vesselboost '; osf -p abk4p fetch /pretrained_models/manual_0429 ' vesselBoostModel],'-echo');
        %%%% ants
        src.ants = 'ml ants/2.5.3';
        system([src.ants '; antsRegistration --version > /dev/null'],'-echo');
        %%%% freesurfer
        src.fs   = 'ml freesurfer/8.0.0';
        system([src.fs   '; mri_convert > /dev/null'],'-echo');
        %%%% vmtk
        src.vmtk = 'ml vmtk/1.5.0';
        system([src.vmtk '; vmtkcenterlines --help > /dev/null'],'-echo');
        %%%% nipype
        src.nipype = 'ml nipype/1.8.3';
        system([src.nipype '; python -c "import nibabel, skimage.morphology; print(''nipype OK'')"'],'-echo');
    case 2
        warning('neurodesk not implemented for this environment');
    otherwise
        dbstack; error('not implemented')
        % neurodeskModule = {
        % ":/neurodesktop-storage/containers/freesurfer_8.0.0_20250210"
        % ":/neurodesktop-storage/containers/afni_24.3.00_20241003"};
        % for i = 1:length(neurodeskModule)
        %     if contains(getenv("PATH"),neurodeskModule{i}); continue; end
        %     setenv("PATH",getenv("PATH") + neurodeskModule{i});
        % end
end

% %%% conda (base system)
% switch envId
%     case 1
%         %%%% skan (https://skeleton-analysis.org/stable/getting_started/install.html)
%         [envNotFound,result] = system('ml nipype/1.8.3; python -c "import skan; print(''skan imported successfully'')"');
%         if envNotFound; disp(['skan not installed. Run:' newline 'conda install -c conda-forge skan' newline 'to install it on your base system.']); else; disp('skan already installed'); end
%     otherwise
%         dbstack; error('not implemented')
% end
    

%% %%%%%%%%%%%%%%%%%%
disp(projectCode)
disp(projectStorage)
disp(projectScratch)
scriptDir = fullfile(projectCode,scriptDir); if ~exist(scriptDir,'dir'); mkdir(scriptDir); end





%%%%%%%%%%%%%%%%%%%
%% Get data pointer
%%%%%%%%%%%%%%%%%%%
% sub-2_vesselPointer20260124161412.mat

fName = 'sub-2_vesselPointer.mat';
% S=2; arxvDir = '/local/users/Proulx-S/db/vsmDiamCenSur/sub-vsmDiamCenSurP2_acq-vfMRI_prsc-dflt_venc-none';
% dir(fullfile(projectCode,'doIt_singleSlabTofVolPrc',['sub-' num2str(S) '_vesselPointer*.mat']))
% load(fullfile(projectCode,'doIt_singleSlabTofVolPrc',['sub-' num2str(S) '_vesselPointer20260125165227.mat']));
load(fullfile(projectCode,'doIt_singleSlabTofVolPrc',fName));
%% %%%%%%%%%%%%%%%%
vessels = info.vessel;
info.project.code

return


% choose good vessel
vIdxList = [6 7 8 11 12 13 14 17 19];
nVessel = length(vIdxList);
[p,~,~] = fileparts(vessels(vIdxList(1)).mask.f);

tmpDir = fullfile(info.project.code,'tmp'); if ~exist(tmpDir,'dir'); mkdir(tmpDir); end

fMaskList = [vessels(vIdxList).mask]; fMaskList = {fMaskList.f}';
fTofList  = [vessels(vIdxList).tof]; fTofList = {fTofList.f}';
mri1 = MRIread(fTofList{1}{1});
mri2 = MRIread(fTofList{1}{2});
mri2.volsize./mri1.volsize

% fTof = info.subject.tof.fList{end};

vmtk_viewVol(fTofList{1}{1});

fSurfList              = cell(nVessel, 1);
fCenterlineList        = cell(nVessel, 1);
fCenterlineRefinedList = cell(nVessel, 1);
for v = 1:nVessel
    vIdx = vIdxList(v);
    fSurfList{v}              = replace(fMaskList{v}, '_mask.nii.gz'   , '_surf.vtk');
    fCenterlineList{v}        = replace(fMaskList{v}, '_mask.nii.gz'   , '_centerline.vtk');
    fCenterlineRefinedList{v} = replace(fMaskList{v}, '_mask.nii.gz'   , '_centerlineRefined.vtk');
end

% Create surface for each vessel: marching cubes -> clean -> smoothing [-> upsample]
forceThis = 1;
for v = 1:nVessel
    if forceThis || ~exist(fSurfList{v},'file')
        vmtk_surfFromSeg(fMaskList{v}, fSurfList{v});
        vmtk_surfClean(fSurfList{v}, fSurfList{v});
        options.passband = 0.01;
        options.iterations = 2000;
        vmtk_surfSmoothing(fSurfList{v}, fSurfList{v}, options);

        vmtk_viewVolAndSurf(fMaskList{v}, fSurfList{v});
        vmtk_viewVolAndSurf(fTofList{v}{1}, fSurfList{v});

    end
end

% Extract centerlines from each surface
for v = 1:nVessel
    if forceThis || ~exist(fCenterlineList{v},'file')
        cmd = vmtk_centerlinesFromSurf(fSurfList{v}, fCenterlineList{v}, [], [], 1);
        vmtk_viewVolAndSurf(fMaskList{v}, fCenterlineList{v});
        vmtk_viewVolAndSurf(fTofList{v}{1}, fCenterlineList{v});
        vmtk_viewVolAndSurf(fTofList{v}{2}, fCenterlineList{v});

        % vmtk_viewVolAndSurf(fTof, fCenterlineList{v});
        % vmtk_viewVolAndSurf(fSeg, fCenterlineList{v});
        % vmtk_viewVolAndSurf(info.vessel(vIdx).seg.fList{end-2}, fCenterlineList{v});
        
    end
end

% Visualize volume with surface or centerlines (uncomment and set v as needed)
% v = 1;
% vmtk_viewVolAndSurf(scaleMaxTof, vesselSurfList{v});
% vmtk_viewVolAndSurf(scaleMaxTof, vesselCenterlineList{v});
% % If vmtk_viewVolAndCenterlines exists in your bassWrap-reg:
% % vmtk_viewVolAndCenterlines(scaleMaxTof, vesselCenterlineList{v}, 1);

% Active tube refinement
opts = struct('iterations', 100, 'potentialweight', 1.0, 'stiffnessweight', 1.0, 'forceThis', forceThis);
for v = 1:nVessel
    if forceThis || ~exist(fCenterlineRefinedList{v},'file')
        fCenterlineRefinedList{v} = vmtk_activetubes(scaleMaxTof, fCenterlineList{v}, fCenterlineRefinedList{v}, opts);
    end
end


