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


if 0
    %%%%%%%%%%%%%%%%%%
    %% Copy data files
    %%%%%%%%%%%%%%%%%
    % get data pointer
    tmp = fullfile(storageDrive,'vsmCenSur','dataPointers.mat');
    pointerFile = fullfile(projectScratch,'dataPointers.mat');
    copyfile(tmp,pointerFile);

    % get nifis
    load(pointerFile);
    s = 1; % perferct all over
    s = 2; % somewhate rigid movement mostly in the first and 2nd run -> mostly correctable with matlab registration over vesselRegion, but frames failed (will probably require temporal smoothing)
    s = 3; % ok
    s = 4; % minimal possibly non-rigid movement
    s = 5; % some non-rigid movement particularly in one vessel on the left on the last run
    s = 6; % ok
    s = 7; % some movements, not clear if rigid

    %tof
    in = fullfile(roi{S}.(acq).(task).rCond.volAnat.tof.folder,roi{S}.(acq).(task).rCond.volAnat.tof.name);
    out = fullfile(projectCode,'data','tof'); if ~exist(out,'dir'); mkdir(out); end
    tof = fullfile(out,'tof.nii.gz');
    copyfile(in,tof);
    switch S
        case 2
            in = fullfile(rxivDir,'tofBrain.nii.gz');
            out = fullfile(projectCode,'data','tofBrain'); if ~exist(out,'dir'); mkdir(out); end
            tofBrain = fullfile(out,'tof.nii.gz');
            copyfile(in,tofBrain);
        otherwise
            error('not implemented')
    end
    

    %vfMRI
    in = roi{S}.(acq).(task).rCond.fOrigList{1};
    out = fullfile(projectCode,'data','vfMRI'); if ~exist(out,'dir'); mkdir(out); end
    vfMRI = fullfile(out,'vfMRI.nii.gz');
    copyfile(in,vfMRI);
    %% %%%%%%%%%%%%%%%
else
    %%%%%%%%%%%%%%%%%
    %% Get data files
    tof      = fullfile(projectCode,'data','tof'     ,'tof.nii.gz'  );
    tofBrain = fullfile(projectCode,'data','tofBrain','tof.nii.gz'  );
    vfMRI    = fullfile(projectCode,'data','vfMRI'   ,'vfMRI.nii.gz');
    %% %%%%%%%%%%%%%%
end



forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vesselboost -- on original resolution tof (with preprocessing, mostly just for the preprocessing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find brain bounding box
mri = MRIread(tofBrain);
[dim1, dim2, dim3] = ind2sub(size(mri.vol), find(mri.vol > 0));
dim1 = [min(dim1) max(dim1)]; if mod(mean(dim1),1); dim1(1) = dim1(1) - 1; end; dim1 = dim1 + [-1 1];
dim2 = [min(dim2) max(dim2)]; if mod(mean(dim2),1); dim2(1) = dim2(1) - 1; end; dim2 = dim2 + [-1 1];
dim3 = [min(dim3) max(dim3)];
cropCenter = [mean(dim2) mean(dim1) 0];
cropSize = [range(dim2)+1 range(dim1)+1 range(dim3)+1];
% Crop to brain bounding box
tofCrop = replace(tof, '/data/tof/tof.nii.gz', ['/tmp/crop/tof.nii.gz']);
if ~exist(fileparts(tofCrop),'dir'); mkdir(fileparts(tofCrop)); end
if forceThis || ~exist(tofCrop,'file')
    cmd = {src.fs};
    cmd{end+1} = 'mri_convert \';
    cmd{end+1} = ['--crop '     strjoin(arrayfun(@num2str, cropCenter, 'UniformOutput', false), ' ') ' \'];
    cmd{end+1} = ['--cropsize ' strjoin(arrayfun(@num2str, cropSize  , 'UniformOutput', false), ' ') ' \'];
    cmd{end+1} = [tof ' \'];
    cmd{end+1} = [tofCrop];
    system(strjoin(cmd,newline),'-echo');
end
% Vesselboost (mostly just for the preprocessing)
tofCrop                                                          ;
tofCropPrc    = replace(tofCrop   , '/crop/'     , '/cropPrc/'  );
tofCropPrcSeg = replace(tofCropPrc, '/tof.nii.gz', '/seg.nii.gz');
if ~exist(fileparts(tofCropPrc),'dir')   ; mkdir(fileparts(tofCropPrc)); end
if ~exist(fileparts(tofCropPrcSeg),'dir'); mkdir(fileparts(tofCropPrcSeg)); end
if forceThis || ~exist(tofCropPrcSeg,'file') || ~exist(tofCropPrc,'file')
    vesselboost_prediction(...
        tofCrop,...
        tofCropPrcSeg,...
        tofCropPrc,...
        vesselBoostModel,3);
end
% % Copy results
% info.subject.seg.fList = {fullfile(info.project.code,'result','seg','scale-1_seg.nii.gz')};
% if ~exist(fileparts(info.subject.seg.fList{1}),'dir'); mkdir(fileparts(info.subject.seg.fList{1})); end
% if forceThis || ~exist(info.subject.seg.fList{1},'file')
%     copyfile(tofCropPrcSeg,info.subject.seg.fList{1});
% end
% info.subject.tof.fList = {fullfile(info.project.code,'result','tof','scale-1_tof.nii.gz')};
% if ~exist(fileparts(info.subject.tof.fList{1}),'dir'); mkdir(fileparts(info.subject.tof.fList{1})); end
% if forceThis || ~exist(info.subject.tof.fList{1},'file')
%     copyfile(tofCropPrc,info.subject.tof.fList{1});
% end    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofCrop;
tofCropPrc;




forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vesselboost -- on preprocessed tof resampled to different scales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scaleList = [1 2 4 8]; % first in the list must be 1
% Get original high-precision voxel spacing
resOrig = MRIread(tofCropPrc,1);
szOrig  = resOrig.volsize;
resOrig = resOrig.volres;
szOrigAnts  = szOrig( [2 1 3]); % row major -> column major
resOrigAnts = resOrig([2 1 3]); % row major -> column major
% round off voxel spacing for consistency across scales
dimStrSpacing = sprintf('%.*gx%.*gx%.*g', 7, resOrigAnts(1), 7, resOrigAnts(2), 7, resOrigAnts(3));
disp(dimStrSpacing)
resAnts = str2double(strsplit(dimStrSpacing,'x'));
res     = resAnts([2 1 3]);

scaleSpacingStr = cell(1,length(scaleList));
tofCropPrcScl    = cell(1,length(scaleList));
tofCropPrcSclSeg = cell(1,length(scaleList));
for i = 1:length(scaleList)

    % Resample tof to different scales
    if i==1
        tofCropPrcScl{i} = fullfile(fileparts(tofCropPrc), ['scale-' num2str(scaleList(i))    '_inter-nn_tof.nii.gz']);
    else
        tofCropPrcScl{i} = fullfile(fileparts(tofCropPrc), ['scale-' num2str(scaleList(i)) '_inter-cubic_tof.nii.gz']);
    end
    if ~exist(fileparts(tofCropPrcScl{i}),'dir'); mkdir(fileparts(tofCropPrcScl{i})); end
    scaleSpacingStr{i} = sprintf('%.*gx%.*gx%.*g', 16, resAnts(1)./scaleList(i), 16, resAnts(2)./scaleList(i), 16, resAnts(3)./scaleList(i));
    if forceThis || ~exist(tofCropPrcScl{i},'file')
        disp(['upsampling tof (scale=' num2str(scaleList(i)) ')... computing']);
        if i==1
            % Nearest neighbor resampling at the rounded voxel spacing for scale-1 (for consistency across scales)
            cmd = {src.ants};
            cmd{end+1} = 'ResampleImage 3 \';
            cmd{end+1} = [tofCropPrc ' \'];
            cmd{end+1} = [tofCropPrcScl{i} ' \'];
            cmd{end+1} = [scaleSpacingStr{i} ' \'];
            cmd{end+1} = ['0 4 3 6']; % [0 = voxel spacing input, 4 = bspline interpolation, 3 = cubic spline, 6 = float output type]
            [status,result] = system(strjoin(cmd, newline), '-echo');
        else
            % Cubic interpolation upsampling at multiples of the rounded original voxel spacing (for consistency across scales)
            cmd = {src.ants};
            cmd{end+1} = 'ResampleImage 3 \';
            cmd{end+1} = [tofCropPrcScl{1} ' \'];
            cmd{end+1} = [tofCropPrcScl{i} ' \'];
            cmd{end+1} = [scaleSpacingStr{i} ' \'];
            cmd{end+1} = ['0 4 3 6']; % [0 = voxel spacing input, 4 = bspline interpolation, 3 = cubic spline, 6 = float output type]
            [status,result] = system(strjoin(cmd, newline), '-echo');
        end
        disp(['upsampling tof (scale=' num2str(scaleList(i)) ')... done']);
    else
        disp(['upsampling tof (scale=' num2str(scaleList(i)) ')... already done']);
    end

    % Vesselboost
    f = strsplit(tofCropPrcScl{i},filesep); f = replace(f{end},'tof','seg');
    tofCropPrcSclSeg{i} = fullfile(fileparts(tofCropPrcScl{i}),f);
    if ~exist(fileparts(tofCropPrcSclSeg{i}),'dir'); mkdir(fileparts(tofCropPrcSclSeg{i})); end
    if forceThis || ~exist(tofCropPrcSclSeg{i},'file')
        disp(['vesselboost (scale=' num2str(scaleList(i)) ')... computing']);
        vesselboost_prediction(...
            tofCropPrcScl{i}   ,...
            tofCropPrcSclSeg{i},...
            [],...
            vesselBoostModel,4);
        disp(['vesselboost (scale=' num2str(scaleList(i)) ')... done']);
    else
        disp(['vesselboost (scale=' num2str(scaleList(i)) ')... already done']);
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofCropPrcScl;
tofCropPrcSclSeg;




forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multiscale consensus segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofCropPrcSclSegUps          = cell(1,length(scaleList));
tofCropPrcSclSegUpsSum       = fullfile(info.project.code,'tmp','seg','scale-C_seg.nii.gz');

% Upsample all segmentations to the highest resolution
for i = 1:length(scaleList)
    tofCropPrcSclSegUps{i} = fullfile(info.project.code,'tmp','seg',['scale-' num2str(scaleList(i)) '_seg.nii.gz']);
    if ~exist(fileparts(tofCropPrcSclSegUps{i}),'dir'); mkdir(fileparts(tofCropPrcSclSegUps{i})); end
    if forceThis || ~exist(tofCropPrcSclSegUps{i},'file')
        disp(['upsampling segmentation (scale=' num2str(scaleList(i)) ')... computing']);
        cmd = {src.ants};
        cmd{end+1} = 'ResampleImage 3 \';
        cmd{end+1} = [tofCropPrcSclSeg{i}    ' \'];
        cmd{end+1} = [tofCropPrcSclSegUps{i} ' \'];
        cmd{end+1} = [scaleSpacingStr{end}   ' \'];
        cmd{end+1} = '0 1 2'; % [0 = spacing, 1 = nearest neighbor interpolation, 2 = unsigned char output]
        [status,result] = system(strjoin(cmd, newline), '-echo');
        disp(['upsampling segmentation (scale=' num2str(scaleList(i)) ')... done']);
    else
        disp(['upsampling segmentation (scale=' num2str(scaleList(i)) ')... already done']);
    end
end
% Compute consensus segmentation (sum of segmentations at each scale)
if ~exist(fileparts(tofCropPrcSclSegUpsSum),'dir'); mkdir(fileparts(tofCropPrcSclSegUpsSum)); end
if forceThis || ~exist(tofCropPrcSclSegUpsSum,'file')
    disp('consensus segmentation (sum of segmentations at each scale)... computing');
    mriSum = MRIread(tofCropPrcSclSegUps{i});
    mriSum.fspec = tofCropPrcSclSegUpsSum;
    for i = 2:length(tofCropPrcSclSegUps)
        mri = MRIread(tofCropPrcSclSegUps{i});
        mriSum.vol = mriSum.vol + mri.vol;
    end
    MRIwrite(mriSum, tofCropPrcSclSegUpsSum,'ushort');
    disp('consensus segmentation (sum of segmentations at each scale)... done');
else
    disp('consensus segmentation (sum of segmentations at each scale)... already done');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




forceThis = 0;
%%%%%%%%%%%%%%%%
%% Store results
%%%%%%%%%%%%%%%%
vesselBoostDir = fullfile(info.project.code,'msVesselBoost'); if ~exist(vesselBoostDir,'dir'); mkdir(vesselBoostDir); end

% tof
info.subject.tof.fList = {};
outDir = fullfile(vesselBoostDir,'tof'); if ~exist(outDir,'dir'); mkdir(outDir); end
info.subject.tof.fList{end+1,1} = fullfile(outDir,'scale-1_tof.nii.gz');
if forceThis || ~exist(info.subject.tof.fList{end},'file')
    copyfile(tofCropPrc,info.subject.tof.fList{end});
end
info.subject.tof.fList{end+1,1} = fullfile(outDir,'scale-8_tof.nii.gz');
if forceThis || ~exist(info.subject.tof.fList{end},'file')
    copyfile(tofCropPrcScl{end},info.subject.tof.fList{end});
end

% multiscale and consensus segmentation
info.subject.seg.fList = {};
outDir = fullfile(vesselBoostDir,'seg'); if ~exist(outDir,'dir'); mkdir(outDir); end
for s = 1:length(scaleList)
    info.subject.seg.fList{end+1,1} = fullfile(outDir,['scale-' num2str(scaleList(s)) '_seg.nii.gz']);
    if forceThis || ~exist(info.subject.seg.fList{end},'file')
        copyfile(tofCropPrcSclSeg{s},info.subject.seg.fList{end});
    end
end
info.subject.seg.fList{end+1,1} = fullfile(outDir,'scale-C_seg.nii.gz');
if forceThis || ~exist(info.subject.seg.fList{end},'file')
    copyfile(tofCropPrcSclSegUpsSum,info.subject.seg.fList{end});
end
%% %%%%%%%%%%%%%



forceThis = 0;
%%%%%%%%%%%%%%%%%%%%
%% Visualize results
%%%%%%%%%%%%%%%%%%%%
cmd = {src.fs};
cmd{end+1} = 'freeview \';
valRange = [0 1]; valRange(valRange==0) = valRange(valRange==0) + 0.01;
cmd{end+1} = ['-v ' info.subject.seg.fList{end} ':colormap=heat:heatscale=' num2str(valRange(1)) ',' num2str(valRange(2)) ' \']; % concensus segmentation
for i = 1:length(info.subject.seg.fList(1:end-1))
    cmd{end+1} = ['-v ' info.subject.seg.fList{i} ':colormap=heat:heatscale=' num2str(valRange(1)) ',' num2str(valRange(2)) ':resample=nearest \']; % individual-scale segmentations
end
cmd{end+1} = ['-v ' info.subject.tof.fList{1} ':resample=nearest']; % original resolution tof
info.subject.visualize.f = fullfile(fileparts(fileparts(info.subject.tof.fList{1})),'msVesselBoost_view.cmd');
if forceThis || ~exist(info.subject.visualize.f,'file')
    fid = fopen(info.subject.visualize.f,'w');
    fprintf(fid,'%s\n',strjoin(cmd,newline));
    fclose(fid);
end
info.subject.visualize.cmd = cmd;
%% %%%%%%%%%%%%%%%%%
disp(strjoin(cmd,newline));




forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%
%% Save vessel pointer
%%%%%%%%%%%%%%%%%%%%%%
% dir(fullfile(info.project.code,mfilename,'*.mat'))
% /scratch/users/Proulx-S/vesselReg3/doIt_singleSlabTofVolPrc/sub-2_vesselPointer20260203164528.mat
% info.pointerFile.vessel = fullfile(info.project.code,mfilename,['sub-' num2str(S) '_vesselPointer' datestr(now,'yyyymmddHHMMSS') '.mat']);
info.pointerFile.msVesselBoost = fullfile(info.project.code,'msVesselBoost',['sub-' num2str(S) '_vesselBoostPointer.mat']);
if ~exist(fileparts(info.pointerFile.msVesselBoost),'dir'); mkdir(fileparts(info.pointerFile.msVesselBoost)); end
if forceThis || ~exist(info.pointerFile.msVesselBoost,'file')
    save(info.pointerFile.msVesselBoost,'info');
    disp(['Vessel pointer saved to:' newline info.pointerFile.msVesselBoost]);
end
%% %%%%%%%%%%%%%%%%%%%
info;


skipThis = 1;
%%%%%%%%%%%%
%% Clear tmp
%%%%%%%%%%%%
tmpFiles = [tofCropPrcSclSegUps'
{tofCropPrcSclSegUpsSum}
tofCropPrcScl'
tofCropPrcSclSeg'
{tofCrop}
{tofCropPrc}];
for i = 1:numel(tmpFiles)
    if ~skipThis && exist(tmpFiles{i}, 'file')
        delete(tmpFiles{i});
    end
end