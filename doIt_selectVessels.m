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
toClean = {};





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get msVesselBoost pointers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msVesselBoost_pointerFile = fullfile(info.project.code,'msVesselBoost',['sub-' num2str(S) '_vesselBoostPointer.mat']);
load(msVesselBoost_pointerFile); % overwrites info, but should be fine
%% %%%%%%%%%%%%%%%%%%%%%%%%%%
info;



forceThis = 0;
consensusThreshValue = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Individualize and label vessels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fMask  = fullfile(info.project.code, 'selectVessels', 'mask' , ['scale-C_thresh-' num2str(consensusThreshValue) '_desc-all_mask.nii.gz' ]);
fLabel = fullfile(info.project.code, 'selectVessels', 'label', ['scale-C_thresh-' num2str(consensusThreshValue) '_desc-all_label.nii.gz']);
if ~exist(fileparts(fMask ),'dir'); mkdir(fileparts(fMask )); end
if ~exist(fileparts(fLabel),'dir'); mkdir(fileparts(fLabel)); end
if forceThis || ~exist(fMask,'file') || ~exist(fLabel,'file')
    % Create vessel mask by thresholding consensus segmentation)
    disp('thresholding consensus segmentation... computing');
    mriSum  = MRIread(info.subject.seg.fList{end});
    mriMask = mriSum;
    mriMask.fspec = fMask;
    mriMask.vol = mriSum.vol >= consensusThreshValue;
    % Split individual vessels (connected components)
    disp('splitting individual vessels... computing');
    % compute connected components
    CC = bwconncomp(mriMask.vol,26);
    [~,b] = sort(cellfun('length',CC.PixelIdxList),'descend');
    CC.PixelIdxList = CC.PixelIdxList(b);
    mriLabel = mriMask;
    mriLabel.fspec = fLabel;
    mriLabel.vol = zeros(size(mriLabel.vol),'uint16');
    for v = 1:length(CC.PixelIdxList)
        mriLabel.vol(CC.PixelIdxList{v}) = uint16(v);
    end
    MRIwrite(mriLabel, fLabel,'ushort');
    MRIwrite(mriMask , fMask ,'ushort');
else
    disp('thresholding and splitting vessels... already done');
end
toClean = cat(1,toClean,{fLabel});
toClean = cat(1,toClean,{fMask});
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% info.subject.mask.fList = {fMask};
info.subject.label.fList = {fLabel};


forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare vfMRI data if available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vfMRI_in = dir(fullfile(info.project.code,'data','vfMRI','vfMRI.nii.gz'));
if ~isempty(vfMRI_in)
    % Average the vfMRI timeseries across time using ants' AverageImages
    vfMRI_in  = fullfile(vfMRI_in.folder, vfMRI_in.name);
    vfMRI_out = fullfile(info.project.code, 'selectVessels', 'vfMRI', 'vfMRI.nii.gz');
    if ~exist(fileparts(vfMRI_out),'dir'); mkdir(fileparts(vfMRI_out)); end
    if forceThis || ~exist(vfMRI_out, 'file')
        cmd = {src.ants};
        cmd{end+1} = ['antsMotionCorr -d 3 \']; % 4D image, average across time
        cmd{end+1} = ['-a ' vfMRI_in ' \'];
        cmd{end+1} = ['-o ' vfMRI_out]; % 0 = average (not sum)
        [status, result] = system(strjoin(cmd, newline), '-echo');
    end
    info.subject.vfMRI.fList = {vfMRI_out};
else
    clear vfMRI_in
    warning(['vfMRI file not found: ' vfMRI_file]);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info.subject.vfMRI;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% View vessel label map in freeview and select vessel labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmd = {src.fs};
cmd{end+1} = 'freeview \';
valRange = [1 1];
cmd{end+1} = ['-v ' fMask  ':colormap=jet:colorscaleq=' num2str(valRange(1)) ',' num2str(valRange(2)) ' \']; % concensus segmentation
cmd{end+1} = ['-v ' fLabel ':colormap=lut \'                       ]; % concensus segmentation
cmd{end+1} = ['-v ' info.subject.tof.fList{1} ':resample=nearest \']; % original resolution tof
cmd{end+1} = ['-v ' info.subject.vfMRI.fList{1}                    ]; % vfMRI

info.subject.visualize(end+1).f = fullfile(fileparts(fileparts(fLabel)),'selectVessels_view.cmd');
if forceThis || ~exist(info.subject.visualize(end).f,'file')
    fid = fopen(info.subject.visualize(end).f,'w');
    fprintf(fid,'%s\n',strjoin(cmd,newline));
    fclose(fid);
end
info.subject.visualize(end).cmd = cmd';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info.subject.visualize(end)
disp(strjoin(replace(info.subject.visualize(end).cmd,'/scratch/users/Proulx-S','/scratch/Proulx-S'),newline));
info.subject.label.selectedVesselIdx = sort([7 9 47 19 27 23 17 15 8 16 100 31 28 12 13 48 26 11 20 36 14]);




forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Relabel selected vessels
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fMask  = fullfile(info.project.code, 'selectVessels', 'mask' , ['scale-C_thresh-' num2str(consensusThreshValue) '_desc-selected_mask.nii.gz' ]);
fLabel = fullfile(info.project.code, 'selectVessels', 'label', ['scale-C_thresh-' num2str(consensusThreshValue) '_desc-selected_label.nii.gz']);
if ~exist(fileparts(fMask ),'dir'); mkdir(fileparts(fMask )); end
if ~exist(fileparts(fLabel),'dir'); mkdir(fileparts(fLabel)); end
if forceThis || ~exist(fLabel,'file')% || ~exist(fMask,'file')
    disp('thresholding consensus segmentation... computing');
    if ~exist('mriLabel','var'); mriLabel = MRIread(info.subject.label.fList{1}); end
    % if ~exist('mriMask' ,'var'); mriMask  = MRIread(info.subject.mask.fList{1} ); end
    % write mask and label files with only selected vessels
    % mriMask.vol( ~ismember(mriMask.vol ,info.subject.label.selectedVesselIdx)) = 0;
    mriLabel.vol(~ismember(mriLabel.vol,info.subject.label.selectedVesselIdx)) = 0;
    for v = 1:length(info.subject.label.selectedVesselIdx)
        % mriMask.vol( mriMask.vol ==info.subject.label.selectedVesselIdx(v)) = 1;
        mriLabel.vol(mriLabel.vol==info.subject.label.selectedVesselIdx(v)) = uint16(v);
    end
    disp('saving...');
    % MRIwrite(mriLabel, fLabel,'ushort');
    MRIwrite(mriMask , fMask ,'ushort');
    disp('thresholding and splitting vessels... done');
else
    disp('thresholding and splitting vessels... already done');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%
% info.subject.mask.fList{end+1}  = fMask;
info.subject.label.fList{end+1} = fLabel;




forceThis = 1;
%%%%%%%%%%%%%%%%%%%%%%
%% Save vessel pointer
%%%%%%%%%%%%%%%%%%%%%%
info.pointerFile.selectVessels = fullfile(info.project.code,'selectVessels',['sub-' num2str(S) '_selectVesselsPointer.mat']);
if ~exist(fileparts(info.pointerFile.selectVessels),'dir'); mkdir(fileparts(info.pointerFile.selectVessels)); end
if forceThis || ~exist(info.pointerFile.selectVessels,'file')
    save(info.pointerFile.selectVessels,'info');
    disp(['Vessel pointer saved to:' newline info.pointerFile.selectVessels]);
end
%% %%%%%%%%%%%%%%%%%%%
info;




skipThis = 1;
%%%%%%%%%%%%
%% Clear tmp
%%%%%%%%%%%%
for i = 1:numel(toClean)
    if ~skipThis && exist(toClean{i}, 'file')
        delete(toClean{i});
    end
end
%% %%%%%%%%%

