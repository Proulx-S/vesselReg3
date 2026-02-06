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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vesselboost -- on original resolution tof (with preprocessing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Vesselboost
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
% Copy results
info.subject.seg.fList = {fullfile(info.project.code,'result','seg','scale-1_seg.nii.gz')};
if ~exist(fileparts(info.subject.seg.fList{1}),'dir'); mkdir(fileparts(info.subject.seg.fList{1})); end
if forceThis || ~exist(info.subject.seg.fList{1},'file')
    copyfile(tofCropPrcSeg,info.subject.seg.fList{1});
end
info.subject.tof.fList = {fullfile(info.project.code,'result','tof','scale-1_tof.nii.gz')};
if ~exist(fileparts(info.subject.tof.fList{1}),'dir'); mkdir(fileparts(info.subject.tof.fList{1})); end
if forceThis || ~exist(info.subject.tof.fList{1},'file')
    copyfile(tofCropPrc,info.subject.tof.fList{1});
end    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofCrop;
tofCropPrc;
tofCropPrcSeg;
info.subject.tof.fList;
info.subject.seg.fList;




forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multiscale Vesselboost -- on upsampled preprocessed tof (different scales)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get original voxel size
cmd = {src.fs};
cmd{end+1} = 'mri_info --res \';
cmd{end+1} = tofCropPrc;
[~,res] = system(strjoin(cmd,newline));
res = str2num(res); res = res(1:3);

usList = [2 4 8];
tofCropPrcScl    = cell(length(usList), 1);
tofCropPrcSclSeg = cell(length(usList), 1);
for i = 1:length(usList)
    tofCropPrcScl{i}    = replace(tofCropPrc, '/cropPrc/', ['/cropPrcScale' num2str(usList(i)) '/']);
    if ~exist(fileparts(tofCropPrcScl{i}),'dir'); mkdir(fileparts(tofCropPrcScl{i})); end
    if forceThis || ~exist(tofCropPrcScl{i},'file')
        % Upsample (cubic interpolation) preprocessed tof
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cmd{end+1} = ['--voxsize ' strjoin(arrayfun(@(x) num2str(x,'%.15g'), res./usList(i), 'UniformOutput', false), ' ') ' \'];
        cmd{end+1} = ['--resample_type cubic \'];
        cmd{end+1} = [tofCropPrc ' \'];
        cmd{end+1} = [tofCropPrcScl{i}];
        system(strjoin(cmd,newline),'-echo');
    else
        disp(['upsampling... already done']);
    end
    % Vesselboost
    tofCropPrcSclSeg{i} = replace(tofCropPrcScl{i},'/tof.nii.gz','/seg.nii.gz');
    if ~exist(fileparts(tofCropPrcSclSeg{i}),'dir'); mkdir(fileparts(tofCropPrcSclSeg{i})); end
    if forceThis || ~exist(tofCropPrcSclSeg{i},'file')
        vesselboost_prediction(...
            tofCropPrcScl{i}   ,...
            tofCropPrcSclSeg{i},...
            [],...
            vesselBoostModel,4);
    else
        disp(['vesselboost... already done']);
    end
    % Copy results
    info.subject.seg.fList{end+1} = replace(info.subject.seg.fList{1},'scale-1',['scale-' num2str(usList(i))]);
    if ~exist(fileparts(info.subject.seg.fList{end}),'dir'); mkdir(fileparts(info.subject.seg.fList{end})); end
    if forceThis || ~exist(info.subject.seg.fList{end},'file')
        copyfile(tofCropPrcSclSeg{i},info.subject.seg.fList{end});
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofCropPrcScl;
tofCropPrcSclSeg;
info.subject.seg.fList;



forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multiscale Vesselboost -- consensus segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,b] = max(usList);
maxRes = res./usList(b);
% Upsample original tof to highest resolution -- for later visualization
upsampledTof = replace(tofCropPrc,'/tof.nii.gz','/tofUpsampled.nii.gz');
if ~exist(fileparts(upsampledTof),'dir'); mkdir(fileparts(upsampledTof)); end
if forceThis || ~exist(upsampledTof,'file')
    cmd = {src.fs};
    cmd{end+1} = 'mri_convert \';
    cmd{end+1} = ['--voxsize ' strjoin(arrayfun(@(x) num2str(x,'%.15g'), maxRes, 'UniformOutput', false), ' ') ' \'];
    cmd{end+1} = ['--resample_type nearest \'];
    cmd{end+1} = [tofCropPrc ' \'];
    cmd{end+1} = [upsampledTof];
    system(strjoin(cmd,newline),'-echo');
    disp(['upsampling tof from scale=1 (original resolution) to the highest resolution... done']);
else
    disp(['upsampling tof from scale=1 (original resolution) to the highest resolution... already done']);
end
% Upsample segementation to highest resolution -- original-resolution-derived segmentation
upsampledSegList = cell(length(tofCropPrcSclSeg)+1,1);
disp(['upsampling segmentation from scale=1 (original resolution) to the highest resolution... computing']);
upsampledSegList{1} = replace(tofCropPrcSeg,'/seg.nii.gz','/segUpsampled.nii.gz');
if ~exist(fileparts(upsampledSegList{1}),'dir'); mkdir(fileparts(upsampledSegList{1})); end
if forceThis || ~exist(upsampledSegList{1},'file')
    cmd = {src.fs};
    cmd{end+1} = 'mri_convert \';
    cmd{end+1} = ['--voxsize ' strjoin(arrayfun(@(x) num2str(x,'%.15g'), maxRes, 'UniformOutput', false), ' ') ' \'];
    cmd{end+1} = ['--resample_type nearest \'];
    cmd{end+1} = ['--out_data_type ' 'uchar' ' \'];
    cmd{end+1} = '--no_scale 1 \'; % don't scale data to 0-255 range
    cmd{end+1} = [tofCropPrcSeg ' \'];
    cmd{end+1} = [upsampledSegList{1}];
    system(strjoin(cmd,newline),'-echo');
    disp(['upsampling segmentation from scale=1 (original resolution) to the highest resolution... done']);
else
    disp(['upsampling segmentation from scale=1 (original resolution) to the highest resolution... already done']);
end
% Upsample segementation to highest resolution -- all-resolution-derived segmentation
disp(['upsampling segmentation from resolutions (scale=' strjoin(arrayfun(@num2str, usList, 'UniformOutput', false), ',') ') to the highest resolution... computing']);
for p = 1:length(tofCropPrcSclSeg)
    disp(['scale=' num2str(usList(p)) '... computing']);
    upsampledSegList{p+1} = replace(tofCropPrcSclSeg{p},'/seg.nii.gz','/segUpsampled.nii.gz');
    if ~exist(fileparts(upsampledSegList{p+1}),'dir'); mkdir(fileparts(upsampledSegList{p+1})); end
    if forceThis || ~exist(upsampledSegList{p+1},'file')
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        if usList(p)~=max(usList)
            cmd{end+1} = ['--voxsize ' strjoin(arrayfun(@(x) num2str(x,'%.15g'), maxRes, 'UniformOutput', false), ' ') ' \'];
            cmd{end+1} = ['--resample_type nearest \'];
        end
        cmd{end+1} = ['--out_data_type ' 'uchar' ' \'];
        cmd{end+1} = '--no_scale 1 \'; % don't scale data to 0-255 range
        cmd{end+1} = [tofCropPrcSclSeg{p} ' \'];
        cmd{end+1} = [upsampledSegList{p+1}];
        system(strjoin(cmd,newline),'-echo');
        disp(['scale=' num2str(usList(p)) '... done']);
    else
        disp(['scale=' num2str(usList(p)) '... already done']);
    end
end
% Compute consensus segmentation (add all segmentations)
consensusSeg = replace(tofCropPrcSclSeg{end},['Scale' num2str(usList(end)) '/'],['ScaleC/']);
if ~exist(fileparts(consensusSeg),'dir'); mkdir(fileparts(consensusSeg)); end
disp(['computing consensus (sum) segmentation... computing']);
if forceThis || ~exist(consensusSeg,'file')
    copyfile(upsampledSegList{1},consensusSeg);
    for p = 2:length(upsampledSegList)
        disp(['adding scale=' num2str(usList(p-1)) ' to consensus segmentation... computing']);
        cmd = {src.fs};
        cmd{end+1} = 'mri_concat \';
        cmd{end+1} = [consensusSeg ' ' upsampledSegList{p} ' \'];
        cmd{end+1} = ['--o ' consensusSeg ' \'];
        cmd{end+1} = '--sum --keep-datatype';
        system(strjoin(cmd,newline),'-echo');
        disp(['adding scale=' num2str(usList(p-1)) ' to consensus segmentation... done']);
    end
    disp(['computing consensus (sum) segmentation... done']);
else
    disp(['computing consensus (sum) segmentation... already done']);
end
% Copy results
info.subject.tof.fList{end+1} = replace(info.subject.tof.fList{1},'scale-1','scale-1toC');
if ~exist(fileparts(info.subject.tof.fList{end}),'dir'); mkdir(fileparts(info.subject.tof.fList{end})); end
if forceThis || ~exist(info.subject.tof.fList{end},'file')
    copyfile(upsampledTof,info.subject.tof.fList{end});
end
info.subject.seg.fList{end+1} = replace(info.subject.seg.fList{1},'scale-1','scale-C');
if ~exist(fileparts(info.subject.seg.fList{end}),'dir'); mkdir(fileparts(info.subject.seg.fList{end})); end
if forceThis || ~exist(info.subject.seg.fList{end},'file')
    copyfile(consensusSeg,info.subject.seg.fList{end});
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
upsampledTof
upsampledSegList;
info.subject.tof.fList;
info.subject.seg.fList;





%%%%%%%%%%%%%%%%%%%%
%% Visualize results
%%%%%%%%%%%%%%%%%%%%
cmd = {src.fs};
cmd{end+1} = 'freeview \';
valRange = [0 1]; valRange(valRange==0) = valRange(valRange==0) + 0.01;
cmd{end+1} = ['-v ' info.subject.seg.fList{end} ':colormap=heat:heatscale=' num2str(valRange(1)) ',' num2str(valRange(2)) ' \']; % concensus segmentation
for i = 1:length(info.subject.seg.fList)-1
    cmd{end+1} = ['-v ' info.subject.seg.fList{i} ':colormap=heat:heatscale=' num2str(valRange(1)) ',' num2str(valRange(2)) ':resample=nearest \']; % individual-scale segmentations
end
cmd{end+1} = ['-v ' info.subject.tof.fList{1} ':resample=nearest']; % original resolution tof
%% %%%%%%%%%%%%%%%%%
disp(strjoin(cmd,newline));





forceThis = 0;
% Note: we might want to consider a more strict criteria (connectivity parameter for bwconncomp) to better separate vessels that are close to one another.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Individualize vessels into a single-vessel label map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
segList   = [upsampledSegList; consensusSeg];
fLabelList = cell(length(segList), 1);
for s = 1:length(segList)
    [p,f,e] = fileparts(segList{s});
    fLabelList{s} = fullfile(p, ['labels_' f e]);
    if forceThis || ~exist(fLabelList{s},'file')
        disp('vessel label map... computing');
        mriLabel = MRIread(segList{s});
        if s == length(segList) % special case for consensus segmentation
            mriLabel.vol = mriLabel.vol>1;
        else
            mriLabel.vol = mriLabel.vol>0;
        end
        CC = bwconncomp(mriLabel.vol,26);
        [~,b] = sort(cellfun('length',CC.PixelIdxList),'descend');
        CC.PixelIdxList = CC.PixelIdxList(b);
        mriLabel.fspec = fLabelList{s};
        mriLabel.vol = zeros(size(mriLabel.vol),'uint16');
        for v = 1:length(CC.PixelIdxList)
            mriLabel.vol(CC.PixelIdxList{v}) = uint16(v);
        end
        MRIwrite(mriLabel, fLabelList{s},'short');
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fLabelList;



forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select vessels and create a reduced single-vessel label map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fTof           = tofCropPrcScl{end};
fLabel         = fLabelList{end}; % use consensus-segmentation-based single-vessel label map
fLabelSelected = replace(fLabel, 'labels.nii.gz', 'labelsSelected.nii.gz');
if forceThis || ~exist(fLabelSelected,'file')
    % visualize vessel labels
    cmd = {src.fs};
    cmd{end+1} = 'freeview \';
    cmd{end+1} = ['-v ' fTof                    ' \'];
    cmd{end+1} = ['-v ' vfMRI  ':resample=nearest \'];
    cmd{end+1} = ['-v ' fLabel ':resample=nearest'  ];
    disp(strjoin(cmd,newline));
    % select vessels
    vesselIdx = sort([7 9 45 19 27 23 66 8 16 15 17 375 32 28 12 13 49 26 11 20 36]); %1->supuriously connected to the sagittal sinus
    % rewrite labels with only selected vessels
    mriLabel = MRIread(fLabel);
    mriLabel.fspec = fLabelSelected;
    mriLabel.vol(~ismember(mriLabel.vol,vesselIdx)) = 0;
    for v = 1:length(vesselIdx)
        mriLabel.vol(mriLabel.vol==vesselIdx(v)) = uint16(v);
    end
    MRIwrite(mriLabel, fLabelSelected,'uchar');
    disp('vessel label map... done');
end
% Copy results
segStr = strsplit(fLabel,'Scale');
if length(segStr) == 2; segStr = strsplit(segStr{2},{'/','_'}); segStr = segStr{1}; else; dbstack; error('Failed to parse segmentation string from label file name'); end
[p,f,e] = fileparts(info.subject.seg.fList); p = char(unique(p));
info.subject.label.fList = {fullfile(replace(p,'/seg','/label'),['scale-' segStr '_label.nii.gz'])};
if ~exist(fileparts(info.subject.label.fList{1}),'dir'); mkdir(fileparts(info.subject.label.fList{1})); end
if forceThis || ~exist(info.subject.label.fList{1},'file')
    copyfile(fLabelSelected,info.subject.label.fList{1});
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fLabelSelected;
info.subject.label.fList;



forceThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Crop out each single vessel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fSegList    = [upsampledSegList; consensusSeg];
% fSegListStr = {'scale-1'}; for s = 1:length(usList); fSegListStr{end+1} = ['scale-' num2str(usList(s))]; end; fSegListStr{end+1} = 'scale-C';
if ~exist('mriLabel','var')
    mriLabel = MRIread(fLabelSelected);
end
nVessel = length(unique(mriLabel.vol(:)))-1;
fVesselTofList   = cell(nVessel, 3               );
fVesselSegList   = cell(nVessel, length(fSegList));
fVesselLabelList = cell(nVessel, length(fSegList));
fVesselMaskList  = cell(nVessel, 1               );
fLikeList        = cell(nVessel, 2);
tmpDir = tempname; if ~exist(tmpDir,'dir'); mkdir(tmpDir); end
for v = 1:nVessel
    disp(' ')
    disp('--------------------------------');
    disp('--------------------------------');
    disp(['Vessel ' sprintf('%02d',v) ' (' num2str(v) '/' num2str(nVessel) ') cropping out...']);
    disp('---------------');
    disp('---------------');
    fVesselTofList{v,1} = replace(upsampledTof   , '/tofUpsampled.nii.gz' , ['/vessel-' sprintf('%02d',v) '_scale-' '1'                  '_interp-nn_tof.nii.gz']);
    fVesselTofList{v,2} = replace(upsampledTof   , '/tofUpsampled.nii.gz' , ['/vessel-' sprintf('%02d',v) '_scale-' num2str(max(usList)) '_interp-nn_tof.nii.gz']);
    fVesselTofList{v,3} = replace(upsampledTof   , '/tofUpsampled.nii.gz' , ['/vessel-' sprintf('%02d',v) '_scale-' num2str(max(usList)) '_interp-cubic_tof.nii.gz']);
    fVesselMaskList{v}  = replace(fLabelList{end}, '/labels_seg.nii.gz'   , ['/vessel-' sprintf('%02d',v) '_scale-C'                              '_mask.nii.gz']);
    for s = 1:length(fSegList)
        fVesselLabelList{v,s} = replace(fLabelList{s}, '/labels_', ['/vessel-' sprintf('%02d',v) '_labels_']);
        fVesselSegList{v,s}   = replace(fSegList{s}  , '/seg'    , ['/vessel-' sprintf('%02d',v) '_seg'    ]);
    end
    s = 1; [p,f,e] = fileparts(fLabelList{s});
    fLikeList{v,1} = fullfile(p, ['vessel-' sprintf('%02d',v) '_scale-1_like.nii.gz']);
    s = length(fLabelList)-1; [p,f,e] = fileparts(fLabelList{s});
    fLikeList{v,2} = fullfile(p, ['vessel-' sprintf('%02d',v) '_scale-' num2str(usList(s-1)) '_like.nii.gz']);
    % Create vessel mask based on label map derived from consensus segmentation
    if forceThis || ~exist(fLikeList{v,1},'file') || ~exist(fLikeList{v,2},'file')
        mriLike       = mriLabel;
        mriLike.vol   = mriLike.vol==v;
        mriLike.fspec = fLikeList{v,2};
    end
    % Create like volume at low resolution
    if forceThis || ~exist(fLikeList{v,1},'file')
        disp('creating like volume at low resolution... computing');
        % by downsampling the mask back to the original resolution
        mriLike1 = MRIread(tofCropPrc,1);
        mriLike1.fspec = fLikeList{v,1};
        mriLike1.vol = imresize3(mriLike.vol, mriLike1.volsize, 'box');
        MRIwrite(mriLike1, fLikeList{v,1}, 'uchar');
        % finding its bounding box
        [dim1, dim2, dim3] = ind2sub(size(mriLike1.vol), find(mriLike1.vol));
        dim1 = [min(dim1) max(dim1)];
        dim2 = [min(dim2) max(dim2)];
        dim3 = [min(dim3) max(dim3)];
        % then cropping it
        cmd = {src.fsl};
        cmd{end+1} = 'fslroi \';
        cmd{end+1} = [fLikeList{v,1} ' \'];
        cmd{end+1} = [fLikeList{v,1} ' \'];
        cmd{end+1} = [num2str(dim2(1)-1) ' ' num2str(dim2(2)-dim2(1)+1) ' \'];
        cmd{end+1} = [num2str(dim1(1)-1) ' ' num2str(dim1(2)-dim1(1)+1) ' \'];
        cmd{end+1} = [num2str(dim3(1)-1) ' ' num2str(dim3(2)-dim3(1)+1)];
        % and pad it
        cmd{end+1} = src.ants;
        cmd{end+1} = 'ImageMath 3 \';
        cmd{end+1} = [fLikeList{v,1} ' \'];
        cmd{end+1} = 'PadImage \';
        cmd{end+1} = [fLikeList{v,1} ' \'];
        cmd{end+1} = '1';
        system(strjoin(cmd,newline),'-echo');
        disp('creating like volume at low resolution... done');
    else
        disp('creating like volume at low resolution... already done');
    end
    % Create like volume at high resolution
    if forceThis || ~exist(fLikeList{v,2},'file')
        disp('creating like volume at high resolution... computing');
        MRIwrite(mriLike, fLikeList{v,2},'uchar');
        cmd{end+1} = src.fs; 
        cmd{end+1} = 'mri_convert \';
        cmd{end+1} = ['--voxsize ' strjoin(arrayfun(@(x) num2str(x,'%.15g'), res./max(usList), 'UniformOutput', false), ' ') ' \'];
        cmd{end+1} = ['--resample_type nearest \'];
        cmd{end+1} = [fLikeList{v,1} ' \'];
        cmd{end+1} = fLikeList{v,2};
        system(strjoin(cmd,newline),'-echo');
        disp('creating like volume at high resolution... done');
    else
        disp('creating like volume at high resolution... already done');
    end
    % Resample other data to the geometry defined by each like volume
    cmd = {src.fs};
    % mask
    if forceThis || ~exist(fVesselMaskList{v},'file')
        disp('resampling mask...');
        cmd{end+1} = 'mri_vol2vol \';
        cmd{end+1} = ['--mov ' fVesselMaskList{v} ' \'];
        cmd{end+1} = ['--targ ' fLikeList{v,2} ' \'];
        cmd{end+1} = ['--o ' fVesselMaskList{v} ' \'];
        cmd{end+1} = '--regheader --nearest --keep-precision';
    else
        disp('resampling mask... already done');
    end
    for s = 1:length(fSegList)
        % seg
        if forceThis || ~exist(fVesselSegList{v,s},'file')
            disp(['resampling seg at scale ' num2str(s) 'of ' num2str(length(fSegList)) '...']);
            cmd{end+1} = 'mri_vol2vol \';
            cmd{end+1} = ['--mov ' fSegList{s} ' \'];
            cmd{end+1} = ['--targ ' fLikeList{v,2} ' \'];
            cmd{end+1} = ['--o ' fVesselSegList{v,s} ' \'];
            cmd{end+1} = '--regheader --nearest --keep-precision';
        else
            disp(['resampling seg at scale ' num2str(s) 'of ' num2str(length(fSegList)) '... already done']);
        end
        % label
        if forceThis || ~exist(fVesselLabelList{v,s},'file')
            disp(['resampling label at scale ' num2str(s) 'of ' num2str(length(fSegList)) '...']);
            cmd{end+1} = 'mri_vol2vol \';
            cmd{end+1} = ['--mov ' fLabelList{s} ' \'];
            cmd{end+1} = ['--targ ' fLikeList{v,2} ' \'];
            cmd{end+1} = ['--o ' fVesselLabelList{v,s} ' \'];
            cmd{end+1} = '--regheader --nearest --keep-precision';
        else
            disp(['resampling label at scale ' num2str(s) 'of ' num2str(length(fSegList)) '... already done']);
        end
    end
    % tof
    if forceThis || ~exist(fVesselTofList{v,1},'file')
        disp(['resampling tof at scale ' num2str(1) 'of ' num2str(length(fSegList)) '...']);
        cmd{end+1} = 'mri_vol2vol \';
        cmd{end+1} = ['--mov ' tofCropPrc ' \'];
        cmd{end+1} = ['--targ ' fLikeList{v,1} ' \'];
        cmd{end+1} = ['--o ' fVesselTofList{v,1} ' \'];
        cmd{end+1} = '--regheader --cubic --keep-precision';
    else
        disp(['resampling tof at scale ' num2str(1) 'of ' num2str(length(fSegList)) '... already done']);
    end
    % tof upsampled nearest
    if forceThis || ~exist(fVesselTofList{v,2},'file')
        disp(['resampling tof at scale ' num2str(length(fSegList)) 'of ' num2str(length(fSegList)) '...']);
        cmd{end+1} = 'mri_vol2vol \';
        cmd{end+1} = ['--mov ' tofCropPrcScl{end} ' \'];
        cmd{end+1} = ['--targ ' fLikeList{v,2} ' \'];
        cmd{end+1} = ['--o ' fVesselTofList{v,2} ' \'];
        cmd{end+1} = '--regheader --nearest --keep-precision';
    else
        disp(['resampling tof at scale ' num2str(length(fSegList)) 'of ' num2str(length(fSegList)) '... already done']);
    end
    % tof upsampled cubic
    if forceThis || ~exist(fVesselTofList{v,3},'file')
        disp(['resampling tof at scale ' num2str(length(fSegList)) 'of ' num2str(length(fSegList)) '...']);
        cmd{end+1} = 'mri_vol2vol \';
        cmd{end+1} = ['--mov ' tofCropPrcScl{end} ' \'];
        cmd{end+1} = ['--targ ' fLikeList{v,2} ' \'];
        cmd{end+1} = ['--o ' fVesselTofList{v,3} ' \'];
        cmd{end+1} = '--regheader --cubic --keep-precision';
    else
        disp(['resampling tof at scale ' num2str(length(fSegList)) 'of ' num2str(length(fSegList)) '... already done']);
    end
    if length(cmd) > 1
        disp(['resampling everything...']);
        system(strjoin(cmd,newline),'-echo');
        disp(['resampling everything... done.']);
    else
        disp(['everything already done.']);
    end
    disp('---------------');
    disp('---------------');
    disp(['Vessel ' sprintf('%02d',v) ' (' num2str(v) '/' num2str(nVessel) ') cropping out... done.']);
    disp('--------------------------------');
    disp('--------------------------------');
    disp(' ');
end
% Copy results
for v = 1:nVessel
    %tof (original resolution)
    info.vessel(v).tof.f = {fullfile(info.project.code,'result','tof',['vessel-' sprintf('%02d',v) '_tof.nii.gz'])};
    if ~exist(fileparts(info.vessel(v).tof.f{1}),'dir'); mkdir(fileparts(info.vessel(v).tof.f{1})); end
    if forceThis || ~exist(info.vessel(v).tof.f{1},'file')
        copyfile(fVesselTofList{v,1},info.vessel(v).tof.f{1});
    end
    %tof (upsampled nearest)
    info.vessel(v).tof.f{end+1,1} = fullfile(info.project.code,'result','tof',['vessel-' sprintf('%02d',v) '_scale-' num2str(max(usList)) '_interp-nn_tof.nii.gz']);
    if ~exist(fileparts(info.vessel(v).tof.f{2}),'dir'); mkdir(fileparts(info.vessel(v).tof.f{2})); end
    if forceThis || ~exist(info.vessel(v).tof.f{2},'file')
        copyfile(fVesselTofList{v},info.vessel(v).tof.f{2});
    end
    %tof (upsampled cubic)
    info.vessel(v).tof.f{end+1,1} = fullfile(info.project.code,'result','tof',['vessel-' sprintf('%02d',v) '_scale-' num2str(max(usList)) '_interp-nn_tof.nii.gz']);
    if ~exist(fileparts(info.vessel(v).tof.f{3}),'dir'); mkdir(fileparts(info.vessel(v).tof.f{3})); end
    if forceThis || ~exist(info.vessel(v).tof.f{3},'file')
        copyfile(fVesselTofList{v},info.vessel(v).tof.f{3});
    end
    %mask
    str = strsplit(fVesselMaskList{v},'Scale'); str = str{2}; str = strsplit(str,{'/','_'}); str = str{1};
    info.vessel(v).mask.f = fullfile(info.project.code,'result','mask',['vessel-' sprintf('%02d',v) '_scale-' str '_mask.nii.gz']);
    if ~exist(fileparts(info.vessel(v).mask.f),'dir'); mkdir(fileparts(info.vessel(v).mask.f)); end
    if forceThis || ~exist(info.vessel(v).mask.f,'file')
        copyfile(fVesselMaskList{v},info.vessel(v).mask.f);
    end
    %seg and label
    info.vessel(v).seg.fList   = cell(length(  fVesselSegList(v,:)),1);
    info.vessel(v).label.fList = cell(length(fVesselLabelList(v,:)),1);
    for s = 1:length(fVesselSegList(v,:))
        if s==1
            str = 'scale-1';
        elseif s==length(fVesselSegList(v,:))
            str = 'scale-C';
        else
            str = ['scale-' num2str(usList(s-1))];
        end
        %seg
        info.vessel(v).seg.fList{s} = fullfile(info.project.code,'result','seg',['vessel-' sprintf('%02d',v) '_' str '_seg.nii.gz']);
        if ~exist(fileparts(info.vessel(v).seg.fList{s}),'dir'); mkdir(fileparts(info.vessel(v).seg.fList{s})); end
        if forceThis || ~exist(info.vessel(v).seg.fList{s},'file')
            copyfile(fVesselSegList{v,s},info.vessel(v).seg.fList{s});
        end        
        %label
        info.vessel(v).label.fList{s} = fullfile(info.project.code,'result','label',['vessel-' sprintf('%02d',v) '_' str '_label.nii.gz']);
        if ~exist(fileparts(info.vessel(v).label.fList{s}),'dir'); mkdir(fileparts(info.vessel(v).label.fList{s})); end
        if forceThis || ~exist(info.vessel(v).label.fList{s},'file')
            copyfile(fVesselLabelList{v,s},info.vessel(v).label.fList{s});
        end        
    end
end
rmdir(tmpDir, 's'); clear tmpDir;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
fVesselTofList;
fVesselLabelList;
fVesselSegList;
fVesselMaskList;
info.vessel.tof;
info.vessel.label;
info.vessel.seg;
info.vessel.mask;



forceThis = 1;
%%%%%%%%%%%%%%%%%%%%%%
%% Save vessel pointer
%%%%%%%%%%%%%%%%%%%%%%
% dir(fullfile(info.project.code,mfilename,'*.mat'))
% /scratch/users/Proulx-S/vesselReg3/doIt_singleSlabTofVolPrc/sub-2_vesselPointer20260203164528.mat
% info.pointerFile.vessel = fullfile(info.project.code,mfilename,['sub-' num2str(S) '_vesselPointer' datestr(now,'yyyymmddHHMMSS') '.mat']);
info.pointerFile.vessel = fullfile(info.project.code,mfilename,['sub-' num2str(S) '_vesselPointer.mat']);
if ~exist(fileparts(info.pointerFile.vessel),'dir'); mkdir(fileparts(info.pointerFile.vessel)); end
if forceThis || ~exist(info.pointerFile.vessel,'file')
    save(info.pointerFile.vessel,'info');
    disp(['Vessel pointer saved to:' newline info.pointerFile.vessel]);
end
%% %%%%%%%%%%%%%%%%%%%
info;


return


tmp = [info.vessel.tof];
tmp = cat(1,tmp.f);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VMTK: surfaces and centerlines from segmented volumes (from vesselReg/vesselReg2)
%  Uses @bassWrap-reg/vmtk_surfFromSeg.m, vmtk_surfClean, vmtk_surfSmoothing,
%  vmtk_surfUpSample, vmtk_centerlinesFromSurf.m, vmtk_viewVolAndSurf.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forceThis = 1;
doSurfUpSample = false;  % set true to upsample for thin/coarse meshes (butterfly, 1 pass)

surfDir       = fullfile(projectCode,'result','surf');       if ~exist(surfDir,'dir'); mkdir(surfDir); end
centerlineDir = fullfile(projectCode,'result','centerlines'); if ~exist(centerlineDir,'dir'); mkdir(centerlineDir); end

vesselSurfList       = cell(nVessel, 1);
vesselCenterlineList = cell(nVessel, 1);
for v = 1:nVessel
    vesselSurfList{v}       = fullfile(surfDir,       ['vessel-' sprintf('%02d',v) '.vtk']);
    vesselCenterlineList{v} = fullfile(centerlineDir, ['vessel-' sprintf('%02d',v) '.vtk']);
end

% Create surface for each vessel: marching cubes -> clean -> smoothing [-> upsample]
for v = 1:nVessel
    if forceThis || ~exist(vesselSurfList{v},'file')
        fSurfRaw    = replace(vesselSurfList{v}, '.vtk', '_raw.vtk');
        fSurfCleaned = replace(vesselSurfList{v}, '.vtk', '_cleaned.vtk');
        vmtk_surfFromSeg(fSegList{v}, fSurfRaw);
        vmtk_surfClean(fSurfRaw, fSurfCleaned);
        options.passband = 0.01;
        options.iterations = 2000;
        vmtk_surfSmoothing(fSurfCleaned, vesselSurfList{v}, options);
        % if doSurfUpSample
        %     fSurfUp = replace(vesselSurfList{v}, '.vtk', '_up.vtk');
        %     vmtk_surfUpSample(vesselSurfList{v}, fSurfUp);
        %     movefile(fSurfUp, vesselSurfList{v});
        % end
        % So that vmtk_viewVolAndSurf(volume_nii, surface) aligns when volume_nii differs from fSegList{v}
        vtk_add_nifti_metadata(vesselSurfList{v}, fSegList{v});
        vmtk_viewVolAndSurf(fSegList{v}, vesselSurfList{v});
    end
end

% Extract centerlines from each surface
for v = 1:nVessel
    if forceThis || ~exist(vesselCenterlineList{v},'file')
        cmd = vmtk_centerlinesFromSurf(vesselSurfList{v}, vesselCenterlineList{v}, [], [], 1);
        disp(strjoin(cmd,newline));
    end
end

% Visualize volume with surface or centerlines (uncomment and set v as needed)
% v = 1;
% vmtk_viewVolAndSurf(scaleMaxTof, vesselSurfList{v});
% vmtk_viewVolAndSurf(scaleMaxTof, vesselCenterlineList{v});
% % If vmtk_viewVolAndCenterlines exists in your bassWrap-reg:
% % vmtk_viewVolAndCenterlines(scaleMaxTof, vesselCenterlineList{v}, 1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VMTK activetubes: refine centerlines using intensity (from @bassWrap-reg/vmtk_activetubes.m)
%  Requires initial centerlines (e.g. from vmtk_centerlinesFromSurf above).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forceThis = 0;
activetubesDir = fullfile(projectCode,'result','centerlines_refined');
if ~exist(activetubesDir,'dir'); mkdir(activetubesDir); end

vesselCenterlineRefinedList = cell(nVessel, 1);
opts = struct('iterations', 100, 'potentialweight', 1.0, 'stiffnessweight', 1.0, 'forceThis', forceThis);
for v = 1:nVessel
    vesselCenterlineRefinedList{v} = fullfile(activetubesDir, ['vessel-' sprintf('%02d',v) '_refined.vtk']);
    if forceThis || ~exist(vesselCenterlineRefinedList{v},'file')
        vesselCenterlineRefinedList{v} = vmtk_activetubes(scaleMaxTof, vesselCenterlineList{v}, vesselCenterlineRefinedList{v}, opts);
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANTs registration (template; bassWrap-reg has no ANTs wrapper)
%% Use src.ants (e.g. 'ml ants/2.5.3') and run antsRegistration / antsApplyTransforms.
%% For 3D affine registration, @bassWrap-reg/afni_3dAllineate.m is an alternative.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example: rigid register moving volume to fixed, then apply transform.
% Uncomment and set fFixed, fMoving, regDir as needed.
%
% fFixed  = scaleMaxTof;   % or vfMRI reference volume
% fMoving = fullfile(projectCode,'data','vfMRI','vfMRI.nii.gz');
% regDir  = fullfile(projectCode,'result','ants_reg');
% if ~exist(regDir,'dir'); mkdir(regDir); end
% prefix  = fullfile(regDir,'ants');
%
% cmd = {src.ants};
% cmd{end+1} = ['antsRegistration -d 3 -f ' fFixed ' -m ' fMoving ' -o ' prefix ' -n linear --transform Rigid[0.1] --metric MI[' fFixed ',' fMoving ',1,32] --convergence [1000x500x250x0,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox'];
% system(strjoin(cmd,newline),'-echo');
%
% fMovingReg = fullfile(regDir,'vfMRI_reg.nii.gz');
% cmd = {src.ants};
% cmd{end+1} = ['antsApplyTransforms -d 3 -i ' fMoving ' -o ' fMovingReg ' -r ' fFixed ' -t ' [prefix '0GenericAffine.mat']];
% system(strjoin(cmd,newline),'-echo');


info.vessel.pointer