clear all
close all
clc
% figure('Menu','none','ToolBar','none');

projectName = 'vesselReg3';
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
        %%%% freesurfer
        src.fsl   = 'ml fsl/6.0.7.16';
        system([src.fsl   '; fslroi > /dev/null'],'-echo');
        %%%% afni
        src.afni = 'ml afni/24.3.00';
        system([src.afni '; 3dinfo > /dev/null'],'-echo');
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
S=2; rxivDir = '/local/users/Proulx-S/db/vsmDiamCenSur/sub-vsmDiamCenSurP2_acq-vfMRI_prsc-dflt_venc-none';
info.project.code    = projectCode;
info.project.storage = projectStorage;
info.project.scratch = projectScratch;
info.subject.rxiv    = rxivDir;
info.subject.idx     = S;
info.toClean = {};



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get selectVessels pointers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cropVessels_pointerFile = fullfile(info.project.code,'cropVessels',['sub-' num2str(S) '_cropVesselsPointer.mat']);
load(cropVessels_pointerFile); % overwrites info, but should be fine
info.toClean = {};
%% %%%%%%%%%%%%%%%%%%%%%%%%%%
info;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract vessel centerlines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nVessel = length(info.vessel);
fMaskList              = cell(nVessel, 1);
fTofList               = cell(nVessel, 2);
fSurfList              = cell(nVessel, 1);
fCenterlineList        = cell(nVessel, 1);
fCenterlineRefinedList = cell(nVessel, 1);
for v = 1:nVessel
    fMaskList{v}              = info.vessel(v).mask.fList{1};
    fTofList(v,:)             = info.vessel(v).tof.fList;
    fSurfList{v}              = replace(info.vessel(v).mask.fList{1}, '_mask.nii.gz'   , '_surf.vtk');
    fCenterlineList{v}        = replace(info.vessel(v).mask.fList{1}, '_mask.nii.gz'   , '_centerline.vtk');
    fCenterlineRefinedList{v} = replace(info.vessel(v).mask.fList{1}, '_mask.nii.gz'   , '_centerlineRefined.vtk');
end

% Create surface for each vessel: marching cubes -> clean -> smoothing
forceThis = 0;
for v = 1:nVessel
    if forceThis || ~exist(fSurfList{v},'file')
        vmtk_surfFromSeg(fMaskList{v}, fSurfList{v});
        vmtk_surfClean(fSurfList{v}, fSurfList{v});
        options.passband = 0.01;
        options.iterations = 2000;
        vmtk_surfSmoothing(fSurfList{v}, fSurfList{v}, options);
        
        vmtk_viewVolAndSurf(fMaskList{v}, fSurfList{v});
        vmtk_viewVolAndSurf(fTofList{v,1}, fSurfList{v});
    end
end

% Extract centerlines from each surface
forceThis = 0;
for v = 1:nVessel
    if forceThis || ~exist(fCenterlineList{v},'file')
        % cmd = vmtk_centerlinesFromSurf(fSurfList{v}, fCenterlineList{v}, [], [], 1);
        cmd = {src.vmtk};
        cmd{end+1} = 'vmtkcenterlines \';
        cmd{end+1} = ['-ifile ' fSurfList{v} ' \'];
        cmd{end+1} = ['-ofile ' fCenterlineList{v}];
        disp(strjoin(cmd,newline));
    end
end
for v = 1:nVessel
    if exist(fCenterlineList{v},'file')
        vmtk_viewVolAndSurf(fTofList{v,1}, fCenterlineList{v});
    end
end

% Active tube refinement
forceThis = 1;
for v = 1:nVessel
    if forceThis || ~exist(fCenterlineRefinedList{v},'file')
        % fCenterlineRefinedList{v} = vmtk_activetubes(fTofList{v}, fCenterlineList{v}, fCenterlineRefinedList{v}, opts);
        if exist(fCenterlineList{v},'file')
            disp([newline 'vmtk_activetubes: ' fCenterlineRefinedList{v}]);
            cmd = {src.vmtk};
            cmd{end+1} = 'vmtkactivetubes \';
            cmd{end+1} = ['-imagefile ' fTofList{v,1} ' \'];
            cmd{end+1} = ['-ifile ' fCenterlineList{v} ' \'];
            cmd{end+1} = ['-ofile ' fCenterlineRefinedList{v}];
            [status, result] = system(strjoin(cmd,newline),'-echo');
            if status ~= 0; error('vmtk_activetubes failed: %s', result); end
        end
    end
end
for v = 1:nVessel
    if exist(fCenterlineRefinedList{v},'file')
        info.vessel(v).ok = true;
        info.vessel(v).centerline.fList{1} = fCenterlineRefinedList{v};
    else
        info.vessel(v).ok = false;
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%

vIdx = [info.vessel.ok];
vessel = info.vessel(vIdx);

fTofList         = [vessel.tof]; fTofList = [fTofList.fList]'; fTofList = fTofList(:,1);
fCenterlineList  = [vessel.centerline]; fCenterlineList = [fCenterlineList.fList]';
nVessel = size(fCenterlineList,1);


fTofRef = info.subject.tof.fList{1};
fCenterlineList
fTofList






