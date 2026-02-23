clear all
close all
% clc
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
selectVessels_pointerFile = fullfile(info.project.code,'selectVessels',['sub-' num2str(S) '_selectVesselsPointer.mat']);
load(selectVessels_pointerFile); % overwrites info, but should be fine
info.toClean = {};
%% %%%%%%%%%%%%%%%%%%%%%%%%%%
info;
nVessel = length(info.subject.label.selectedVesselIdx);



forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write temporary crop masks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fCropMask1 = cell(nVessel, 1);
fCropMaskC = cell(nVessel, 1);
for v = 1:nVessel
    fCropMask1{v} = fullfile(info.project.code, 'tmp', 'cropMask', ['vessel-' sprintf('%02d',v) '_scale-1_cropMask.nii.gz']);
    fCropMaskC{v} = fullfile(info.project.code, 'tmp', 'cropMask', ['vessel-' sprintf('%02d',v) '_scale-C_cropMask.nii.gz']);    
    if forceThis || ~exist(fCropMask1{v}, 'file') || ~exist(fCropMaskC{v}, 'file')
        disp(['writing temporary crop masks for vessel (' num2str(v) '/' num2str(nVessel) ')... computing']);
        if ~exist('mriLabel','var'); mriLabel = MRIread(info.subject.label.fList{2}); end
        if ~exist('mriCrop1','var'); mriCropC = mriLabel; end
        if ~exist('mriCrop1','var'); mriCrop1 = MRIread(info.subject.tof.fList{1},1); end

        mriCrop1.fspec = fCropMask1{v};
        mriCropC.fspec = fCropMaskC{v};
        mriCropC.vol = mriLabel.vol==v;
        mriCrop1.vol = imresize3(mriCropC.vol, mriCrop1.volsize, 'box'    );
        mriCropC.vol = imresize3(mriCrop1.vol, mriCropC.volsize, 'nearest');

        if ~exist(fileparts(fCropMask1{v}), 'dir'); mkdir(fileparts(fCropMask1{v})); end
        if ~exist(fileparts(fCropMaskC{v}), 'dir'); mkdir(fileparts(fCropMaskC{v})); end
        MRIwrite(mriCrop1, fCropMask1{v}, 'ushort');
        MRIwrite(mriCropC, fCropMaskC{v}, 'ushort');
    else
        disp(['writing temporary crop masks for vessel (' num2str(v) '/' num2str(nVessel) ')... already done']);
    end
    info.toClean{end+1,1} = fCropMask1{v};
    info.toClean{end+1,1} = fCropMaskC{v};
end

[~,scaleList,~] = fileparts(replace(replace(info.subject.seg.fList,'.nii',''),'.gz','')); for s = 1:length(info.subject.seg.fList); scaleList{s} = strsplit(scaleList{s}, '_'); scaleList{s} = scaleList{s}{contains(scaleList{s}, 'scale-')}; scaleList{s} = strsplit(scaleList{s},'-'); scaleList{s} = scaleList{s}{end}; end
   scaleMax     = max(cellfun(@str2num, scaleList(1:end-1)));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%
fCropMask1;
fCropMaskC;


forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Crop out each single vessel from tof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fTofList = info.subject.tof.fList;
fVesselTofList = cell(nVessel, length(fTofList));
for v = 1:nVessel
    for s = 1:length(fTofList)
        f = strsplit(fTofList{s}, filesep);
        fVesselTofList{v,s} = fullfile(info.project.code, 'cropVessels', 'tof', ['vessel-' sprintf('%02d',v) '_' f{end}]);
        if ~exist(fileparts(fVesselTofList{v,s}), 'dir'); mkdir(fileparts(fVesselTofList{v,s})); end
        if forceThis || ~exist(fVesselTofList{v,s}, 'file')
            disp(['cropping vessel (' num2str(v) '/' num2str(nVessel) ') from tof (' num2str(s) '/' num2str(length(fTofList)) ')... computing']);
            cmd = {src.ants};
            cmd{end+1} = 'ExtractRegionFromImageByMask 3 \';
            cmd{end+1} = [fTofList{s} ' \'];
            cmd{end+1} = [fVesselTofList{v,s} ' \'];
            scaleStr = strsplit(fTofList{s}, filesep); scaleStr = scaleStr{end}; scaleStr = strsplit(scaleStr, '_'); scaleStr = scaleStr{contains(scaleStr, 'scale-')};
            switch scaleStr
                case 'scale-1'
                    cmd{end+1} = [fCropMask1{v} ' \'];
                    cmd{end+1} = ['1 1']; % [label, padRadius]
                case {['scale-' num2str(scaleMax)] 'scale-C'}
                    cmd{end+1} = [fCropMaskC{v} ' \'];
                    cmd{end+1} = ['1 ' num2str(scaleMax)]; % [label, padRadius]
                otherwise
                    error('Invalid scale');
            end
            [status, result] = system(strjoin(cmd, newline)); if status ~= 0; dbstack; error(result); end
            disp(['cropping vessel (' num2str(v) '/' num2str(nVessel) ') from tof (' num2str(s) '/' num2str(length(fTofList)) ')... done']);
        else
            disp(['cropping vessel (' num2str(v) '/' num2str(nVessel) ') from tof (' num2str(s) '/' num2str(length(fTofList)) ')... already done']);
        end
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fVesselTofList;




forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Crop out each single vessel from seg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scaleInd     = [1 length(scaleList)+[-1 0]];
fSegList     = info.subject.seg.fList(scaleInd);
curScaleList = scaleList(scaleInd);
fVesselSegList = cell(nVessel, length(curScaleList));
for v = 1:nVessel
    for s = 1:length(curScaleList)
        f = strsplit(fSegList{s}, filesep);
        fVesselSegList{v,s} = fullfile(info.project.code, 'cropVessels', 'seg', ['vessel-' sprintf('%02d',v) '_' f{end}]);
        if ~exist(fileparts(fVesselSegList{v,s}), 'dir'); mkdir(fileparts(fVesselSegList{v,s})); end
        if forceThis || ~exist(fVesselSegList{v,s}, 'file')
            disp(['cropping vessel (' num2str(v) '/' num2str(nVessel) ') from seg (' num2str(s) '/' num2str(length(fSegList)) ')... computing']);
            cmd = {src.ants};
            cmd{end+1} = 'ExtractRegionFromImageByMask 3 \';
            cmd{end+1} = [fSegList{s} ' \'];
            cmd{end+1} = [fVesselSegList{v,s} ' \'];
            scaleStr = strsplit(fSegList{s}, filesep); scaleStr = scaleStr{end}; scaleStr = strsplit(scaleStr, '_'); scaleStr = scaleStr{contains(scaleStr, 'scale-')};
            switch scaleStr
                case 'scale-1'
                    cmd{end+1} = [fCropMask1{v} ' \'];
                    cmd{end+1} = ['1 1']; % [label, padRadius]
                case {['scale-' num2str(scaleMax)] 'scale-C'}
                    cmd{end+1} = [fCropMaskC{v} ' \'];
                    cmd{end+1} = ['1 ' num2str(scaleMax)]; % [label, padRadius]
                otherwise
                    error('Invalid scale');
            end
            [status, result] = system(strjoin(cmd, newline)); if status ~= 0; dbstack; error(result); end
            disp(['cropping vessel (' num2str(v) '/' num2str(nVessel) ') from seg (' num2str(s) '/' num2str(length(fSegList)) ')... done']);
        else
            disp(['cropping vessel (' num2str(v) '/' num2str(nVessel) ') from seg (' num2str(s) '/' num2str(length(fSegList)) ')... already done']);
        end
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fVesselSegList;



forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Crop out each single vessel from label
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fLabel = info.subject.label.fList{2};
fVesselLabelList = cell(nVessel, 1);
for v = 1:nVessel
    f = strsplit(fLabel, filesep);
    fVesselLabelList{v} = fullfile(info.project.code, 'cropVessels', 'label', ['vessel-' sprintf('%02d',v) '_' f{end}]);
    if ~exist(fileparts(fVesselLabelList{v}), 'dir'); mkdir(fileparts(fVesselLabelList{v})); end
    if forceThis || ~exist(fVesselLabelList{v}, 'file')
        disp(['cropping vessel (' num2str(v) '/' num2str(nVessel) ') from label... computing']);
        cmd = {src.ants};
        cmd{end+1} = 'ExtractRegionFromImageByMask 3 \';
        cmd{end+1} = [fLabel ' \'];
        cmd{end+1} = [fVesselLabelList{v} ' \'];
        cmd{end+1} = [fCropMaskC{v} ' \'];
        cmd{end+1} = ['1 ' num2str(scaleMax)]; % [label, padRadius]
        [status, result] = system(strjoin(cmd, newline)); if status ~= 0; dbstack; error(result); end
        disp(['cropping vessel (' num2str(v) '/' num2str(nVessel) ') from label... done']);
    else
        disp(['cropping vessel (' num2str(v) '/' num2str(nVessel) ') from label... already done']);
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fVesselLabelList;



forceThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate mask from each vessel label
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fLabel = info.subject.label.fList{2};
fVesselMaskUncroppedList = cell(nVessel, 1);
for v = 1:nVessel
    f = strsplit(fVesselLabelList{v}, filesep);
    f = replace(f{end}, '_label.', '_mask.');
    f = replace(f, 'desc-selected_', 'desc-selectedNotCropped_');
    fVesselMaskUncroppedList{v} = fullfile(info.project.code, 'cropVessels', 'mask', f);
    if ~exist(fileparts(fVesselMaskUncroppedList{v}), 'dir'); mkdir(fileparts(fVesselMaskUncroppedList{v})); end
    if forceThis || ~exist(fVesselMaskUncroppedList{v}, 'file')
        disp(['generating mask from uncropped label for vessel (' num2str(v) '/' num2str(nVessel) ')... computing']);
        cmd = {src.ants};
        cmd{end+1} = ['ThresholdImage 3 \'];
        cmd{end+1} = [fLabel ' \'];
        cmd{end+1} = [fVesselMaskUncroppedList{v} ' \'];
        cmd{end+1} = [num2str(v) ' ' num2str(v) ' 1 0']; % threshlo threshhi insideValue outsideValue
        [status, result] = system(strjoin(cmd, newline)); if status ~= 0; dbstack; error(result); end
        disp(['generating mask from uncropped label for vessel (' num2str(v) '/' num2str(nVessel) ')... done']);
    else
        disp(['generating mask from uncropped label for vessel (' num2str(v) '/' num2str(nVessel) ')... already done']);
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fVesselMaskUncroppedList;




forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate mask from each croppped label
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fVesselMaskList = cell(size(fVesselLabelList));
for v = 1:nVessel
    f = strsplit(fVesselLabelList{v}, filesep);
    f = replace(f{end}, '_label.', '_mask.');
    fVesselMaskList{v} = fullfile(info.project.code, 'cropVessels', 'mask', f);
    if ~exist(fileparts(fVesselMaskList{v}), 'dir'); mkdir(fileparts(fVesselMaskList{v})); end
    if forceThis || ~exist(fVesselMaskList{v}, 'file')
        disp(['generating mask from label for vessel (' num2str(v) '/' num2str(nVessel) ')... computing']);
        cmd = {src.ants};
        cmd{end+1} = ['ThresholdImage 3 \'];
        cmd{end+1} = [fVesselLabelList{v} ' \'];
        cmd{end+1} = [fVesselMaskList{v} ' \'];
        cmd{end+1} = [num2str(v) ' ' num2str(v) ' 1 0']; % threshlo threshhi insideValue outsideValue
        [status, result] = system(strjoin(cmd, newline)); if status ~= 0; dbstack; error(result); end
        disp(['generating mask from label for vessel (' num2str(v) '/' num2str(nVessel) ')... done']);
    else
        disp(['generating mask from label for vessel (' num2str(v) '/' num2str(nVessel) ')... already done']);
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fVesselMaskList;



forceThis = 1;
%%%%%%%%%%%%%%%%%%%%%%
%% Save vessel pointer
%%%%%%%%%%%%%%%%%%%%%%
for v = 1:nVessel
    info.vessel(v,1).tof.fList   = fVesselTofList(v,:)';
    info.vessel(v,1).seg.fList   = fVesselSegList(v,:)';
    info.vessel(v,1).label.fList = fVesselLabelList(v,:)';
    info.vessel(v,1).mask.fList  = fVesselMaskList(v,:)';
    info.vessel(v,1).maskUncropped.fList = fVesselMaskUncroppedList(v,:)';
end
info.pointerFile.cropVessels = fullfile(info.project.code,'cropVessels',['sub-' num2str(S) '_cropVesselsPointer.mat']);
if ~exist(fileparts(info.pointerFile.cropVessels),'dir'); mkdir(fileparts(info.pointerFile.cropVessels)); end
if forceThis || ~exist(info.pointerFile.cropVessels,'file')
    save(info.pointerFile.cropVessels,'info');
    disp(['Crop vessels pointer saved to:' newline info.pointerFile.cropVessels]);
end
%% %%%%%%%%%%%%%%%%%%%
info;




skipThis = 1;
%%%%%%%%%%%%
%% Clear tmp
%%%%%%%%%%%%
for i = 1:numel(info.toClean)
    if ~skipThis && exist(info.toClean{i}, 'file')
        delete(info.toClean{i});
    end
end
%% %%%%%%%%%







fVesselTofList(1,:)'
fVesselSegList(1,:)'
fVesselLabelList(1,:)';
fVesselMaskList(1,:)';
