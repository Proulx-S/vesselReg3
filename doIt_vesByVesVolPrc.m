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
S=2; arxvDir = '/local/users/Proulx-S/db/vsmDiamCenSur/sub-vsmDiamCenSurP2_acq-vfMRI_prsc-dflt_venc-none';
dir(fullfile(projectCode,'doIt_singleSlabTofVolPrc',['sub-' num2str(S) '_vesselPointer*.mat']))
load(fullfile(projectCode,'doIt_singleSlabTofVolPrc',['sub-' num2str(S) '_vesselPointer20260125165227.mat']));
%% %%%%%%%%%%%%%%%%
vessels = info.vessel.pointer;




forceThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute cross-scale consensus segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,prcLabel,~] = fileparts(fileparts(vessels(1).fSeg))
prcInd = contains(prcLabel',{'tofPrcUps' 'tofUps1xPrc'});

for v = 1:length(vessels)
    vLabel = vessels(v).label;
    fSegList          = vessels(v).fSeg(prcInd)';
    fSegList_usFactor = vessels(v).upsampleFactor(prcInd)';
    usRefList          = vessels(v).cropRef.fLabelList';
    usRefList_usFactor = vessels(v).cropRef.usList';

    % Upsample segmentation from all scales to the highest resolution
    [~,b] = max(usRefList_usFactor);
    fLike = usRefList{b};
    fSegListCleaned = cell(length(fSegList),1);
    fSegListCleanedResampled = cell(length(fSegList),1);
    for p = 1:length(fSegList)
        
        % Exclude voxels from other vessels or not connected to the main vessel
        fSegListCleaned{p} = fullfile(scriptDir,['vessel-' sprintf('%02d',vLabel) '_scale-' num2str(fSegList_usFactor(p)) '_mask.nii.gz']);
        if forceThis || ~exist(fSegListCleaned{p},'file')
            mriRef = MRIread(usRefList{fSegList_usFactor(p)==usRefList_usFactor});
            mri    = MRIread(fSegList{p});
            CCRef = bwconncomp(mriRef.vol==vLabel,26);
            CC    = bwconncomp(mri.vol,26);
            % find best match component between the mask upsampled from the vessel of interest and the segmentation on upsampled data
            n = zeros(length(CC.PixelIdxList),1);
            for c = 1:length(CC.PixelIdxList)
                n(c) = nnz(ismember(CCRef.PixelIdxList{1},CC.PixelIdxList{c}));
            end
            [~,b] = max(n);
            % exclude voxels based on that component
            eMask = true(size(mri.vol));
            eMask(CC.PixelIdxList{b}) = false;
            mri.vol(eMask) = 0;
            MRIwrite(mri, fSegListCleaned{p});
        end

        % Upsample the cleaned mask to the highest resolution
        fSegListCleanedResampled{p} = fullfile(scriptDir,['vessel-' sprintf('%02d',vLabel) '_scale-' num2str(fSegList_usFactor(p)) '_maskUpsampled.nii.gz']);
        if forceThis || ~exist(fSegListCleanedResampled{p},'file')
            cmd = {src.fs};
            cmd{end+1} = 'mri_convert \';
            cmd{end+1} = ['--resample_type nearest \'];
            cmd{end+1} = ['--like ' fLike ' \'];
            cmd{end+1} = [fSegListCleaned{p} ' \'];
            cmd{end+1} = [fSegListCleanedResampled{p}];
            system(strjoin(cmd,newline),'-echo');
        end
    end

    % Compute consensus segmentation across scales
    vessels(v).fSegConsensus = fullfile(scriptDir,['vessel-' sprintf('%02d',vLabel) '_scale-consensus.nii.gz']);
    if forceThis || ~exist(vessels(v).fSegConsensus,'file')
        vessel = [];
        for p = 1:length(fSegListCleanedResampled)
            mri = MRIread(fSegListCleanedResampled{p});
            vessel = cat(5,vessel,mri.vol);
        end
        mri.vol = mean(vessel,5);
        MRIwrite(mri, vessels(v).fSegConsensus);
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(vessels(v).fSegConsensus)



