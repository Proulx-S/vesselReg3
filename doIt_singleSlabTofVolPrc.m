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
rxivThis  = 0;
usList = [1 2 3 4 6 8];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vesselboost on cropped -> upsampled -> bias-field-corrected -> denoised tof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofList1 = cell(length(usList), 1);
tofList2 = cell(length(usList), 1);
tofSegList1 = cell(length(usList), 1);
tofSegList2 = cell(length(usList), 1);
for i = 1:length(usList)
    % Crop and upsample tof
    I1 = usList(i);    
    disp('-----------------------');
    disp('-----------------------');
    disp(['Upsampling tof by ' num2str(I1) 'x...']);
    disp('-----------------------');
    disp('-----------------------');
    tof_upsampled = replace(tof, '/data/tof/', ['/tmp/tofUps' num2str(I1) 'x/']);
    tofList1{i} = tof_upsampled;
    if ~exist(fileparts(tof_upsampled),'dir'); mkdir(fileparts(tof_upsampled)); end
    mri = MRIread(tofBrain);
    [dim1, dim2, dim3] = ind2sub(size(mri.vol), find(mri.vol > 0));
    dim1 = [min(dim1) max(dim1)]; if mod(mean(dim1),1); dim1(1) = dim1(1) - 1; end; dim1 = dim1 + [-1 1];
    dim2 = [min(dim2) max(dim2)]; if mod(mean(dim2),1); dim2(1) = dim2(1) - 1; end; dim2 = dim2 + [-1 1];
    dim3 = [min(dim3) max(dim3)];
    if forceThis || ~exist(tof_upsampled,'file')
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cropCenter = [mean(dim2) mean(dim1) 0].*I1;
        cropSize = [range(dim2)+1 range(dim1)+1 range(dim3)+1].*I1;    
        cmd{end+1} = ['--crop '     strjoin(arrayfun(@num2str, cropCenter, 'UniformOutput', false), ' ') ' \'];
        cmd{end+1} = ['--cropsize ' strjoin(arrayfun(@num2str, cropSize  , 'UniformOutput', false), ' ') ' \'];
        if I1~=1
            cmd{end+1} = ['--resample_type cubic \'];
            cmd{end+1} = ['--upsample ' num2str(I1) ' \'];
        end
        cmd{end+1} = [tof ' \'];
        cmd{end+1} = [tof_upsampled];
        system(strjoin(cmd,newline),'-echo');
    end

    % Vesselboost (no preprocessing)
    disp('----------------------------------');
    disp('----------------------------------');
    disp(['Vesselboost (without preprocessing) on upsampled ' num2str(I1) 'x tof...']);
    disp('----------------------------------');
    disp('----------------------------------');

    tof_upsampled_prc = replace(tof_upsampled    , ['/tofUps' num2str(I1) 'x/'], ['/tofUps' num2str(I1) 'x/'   ]); if ~exist(fileparts(tof_upsampled_prc),'dir'); mkdir(fileparts(tof_upsampled_prc)); end
    tof_upsampled_seg = replace(tof_upsampled    , ['/tofUps' num2str(I1) 'x/'], ['/tofUps' num2str(I1) 'xSeg/']); if ~exist(fileparts(tof_upsampled_seg),'dir'); mkdir(fileparts(tof_upsampled_seg)); end
    tofSegList1{i}    = replace(tof_upsampled_seg, '.nii.gz'                   , '.nii'                         );
    if forceThis || ~exist(tofSegList1{i},'file')
        % Run vesselboost
        vesselboost_prediction(fileparts(tof_upsampled),fileparts(tof_upsampled_seg),fileparts(tof_upsampled_prc),vesselBoostModel,4);
        % Reduce precision and decompress
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cmd{end+1} = '--out_data_type uchar \';
        cmd{end+1} = [tof_upsampled_seg ' \'];
        cmd{end+1} = tofSegList1{i};
        system(strjoin(cmd,newline),'-echo');
    end
    % Archive
    rxivFile = strsplit(tofSegList1{i}, '/'); rxivFile = fullfile(rxivDir,[rxivFile{end-1} '.nii.gz']);
    if rxivThis
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cmd{end+1} = [tofSegList1{i} ' \'];
        cmd{end+1} = rxivFile;
        system(strjoin(cmd,newline),'-echo');
    end

    % Vesselboost (including bias-field-correction and denoising preprocessing)
    disp('----------------------------------');
    disp('----------------------------------');
    disp(['Preprocessing and vesselboost on upsampled ' num2str(I1) 'x tof...']);
    disp('----------------------------------');
    disp('----------------------------------');

    tof_upsampled_prc = replace(tof_upsampled    , ['/tofUps' num2str(I1) 'x/'], ['/tofUps' num2str(I1) 'xPrc/'   ]); if ~exist(fileparts(tof_upsampled_prc),'dir'); mkdir(fileparts(tof_upsampled_prc)); end
    tof_upsampled_seg = replace(tof_upsampled    , ['/tofUps' num2str(I1) 'x/'], ['/tofUps' num2str(I1) 'xPrcSeg/']); if ~exist(fileparts(tof_upsampled_seg),'dir'); mkdir(fileparts(tof_upsampled_seg)); end
    tofSegList2{i}    = replace(tof_upsampled_seg, '.nii.gz'                   , '.nii'                            );
    tofList2{i}       = tof_upsampled_prc;
    if forceThis || ~exist(tofSegList2{i},'file')
        % Run vesselboost
        vesselboost_prediction(fileparts(tof_upsampled),fileparts(tof_upsampled_seg),fileparts(tof_upsampled_prc),vesselBoostModel,3);
        % Reduce precision and decompress
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cmd{end+1} = '--out_data_type uchar \';
        cmd{end+1} = [tof_upsampled_seg ' \'];
        cmd{end+1} = tofSegList2{i};
        system(strjoin(cmd,newline),'-echo');
    end
    % Archive
    rxivFile = strsplit(tofSegList2{i}, '/'); rxivFile = fullfile(rxivDir,[rxivFile{end-1} '.nii.gz']);
    if rxivThis
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cmd{end+1} = [tofSegList2{i} ' \'];
        cmd{end+1} = rxivFile;
        system(strjoin(cmd,newline),'-echo');
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(tofList1)
disp(tofList2)
disp(tofSegList1)
disp(tofSegList2)


forceThis = 0;
rxivThis  = 0;
usList = [1 2 3 4 6 8];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vesselboost on cropped -> bias-field-corrected -> denoised -> upsampled tof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofList3    = cell(length(usList), 1);
tofSegList3 = cell(length(usList), 1);
for i = 1:length(usList)
    % Upsample previously cropped, bias-field-corrected and denoised tof
    I1 = usList(i); 
    if I1==1; continue; end

    disp('-----------------------');
    disp('-----------------------');
    disp(['Upsampling preprocessed tof by ' num2str(I1) 'x...']);
    disp('-----------------------');
    disp('-----------------------');
    tof_prc       = replace(tof, '/data/tof/',  '/tmp/tofUps1xPrc/'               );
    tof_upsampled = replace(tof, '/data/tof/', ['/tmp/tofPrcUps' num2str(I1) 'x/']);
    tofList3{i} = tof_upsampled;
    if ~exist(fileparts(tof_upsampled),'dir'); mkdir(fileparts(tof_upsampled)); end
    if forceThis || ~exist(tof_upsampled,'file')
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cmd{end+1} = '--resample_type cubic \';
        cmd{end+1} = ['--upsample ' num2str(I1) ' \'];
        cmd{end+1} = [tof_prc ' \'];
        cmd{end+1} = [tof_upsampled];
        system(strjoin(cmd,newline),'-echo');
    end
    
    % Vesselboost
    disp('----------------------------------');
    disp('----------------------------------');
    disp(['Vesselboost on preprocessed then upsampled ' num2str(I1) 'x tof...']);
    disp('----------------------------------');
    disp('----------------------------------');

    tof_upsampled_prc = replace(tof_upsampled, ['/tofPrcUps' num2str(I1) 'x/'], ['/tofPrcUps' num2str(I1) 'xDum/']); if ~exist(fileparts(tof_upsampled_prc),'dir'); mkdir(fileparts(tof_upsampled_prc)); end
    tof_upsampled_seg = replace(tof_upsampled, ['/tofPrcUps' num2str(I1) 'x/'], ['/tofPrcUps' num2str(I1) 'xSeg/']); if ~exist(fileparts(tof_upsampled_seg),'dir'); mkdir(fileparts(tof_upsampled_seg)); end
    tofSegList3{i}    = replace(tof_upsampled_seg, '.nii.gz', '.nii');
    if forceThis || ~exist(tofSegList3{i},'file')
        % Run vesselboost
        vesselboost_prediction(fileparts(tof_upsampled),fileparts(tof_upsampled_seg),fileparts(tof_upsampled_prc),vesselBoostModel,4);
        % Reduce precision and decompress
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cmd{end+1} = '--out_data_type uchar \';
        cmd{end+1} = [tof_upsampled_seg ' \'];
        cmd{end+1} = tofSegList3{i};
        system(strjoin(cmd,newline),'-echo');
    end

    % Archive
    rxivFile = strsplit(tofSegList3{i}, '/'); rxivFile = fullfile(rxivDir,[rxivFile{end-1} '.nii.gz']);
    if rxivThis
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cmd{end+1} = [tofSegList3{i} ' \'];
        cmd{end+1} = rxivFile;
        system(strjoin(cmd,newline),'-echo');
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(tofList3)
disp(tofSegList3)


tofList = cat(1,...
tofList1(~cellfun('isempty',tofList1)),...
tofList2(~cellfun('isempty',tofList2)),...
tofList3(~cellfun('isempty',tofList3)));
tofSegList = cat(1,...
tofSegList1(~cellfun('isempty',tofSegList1)),...
tofSegList2(~cellfun('isempty',tofSegList2)),...
tofSegList3(~cellfun('isempty',tofSegList3)));



forceThis = 0;
rxivThis  = 0;
%%%%%%%%%%%%%%%%%%%%%%%%
%% Individualize vessels
%%%%%%%%%%%%%%%%%%%%%%%%
% find connected components (vessels) based on segmentation from non-upsampled but preprocessed tof
mriSeg     = MRIread(tofSegList2{1});
mriSeg.vol = logical(mriSeg.vol);
mriTof     = MRIread(replace(tofSegList2{1},'Seg/tof.nii','/tof.nii.gz'));
mriLabels  = mriSeg;
mriLabels.fspec = replace(mriSeg.fspec,'/tof.nii','/tofLabels1.nii');
if forceThis || ~exist(mriLabels.fspec,'file')
    CC = bwconncomp(mriSeg.vol,26);
    % sort vessels from big to small vessel network volume
    [~,b] = sort(cellfun('length',CC.PixelIdxList),'descend');
    CC.PixelIdxList = CC.PixelIdxList(b);
    % label indivudal vessels in a single map for visualization
    mriLabels.vol = zeros(size(mriSeg.vol));
    for v = 1:length(CC.PixelIdxList)
        mriLabels.vol(CC.PixelIdxList{v}) = v;
    end
    MRIwrite(mriLabels, mriLabels.fspec);
else
    mriLabels = MRIread(mriLabels.fspec);
end
% visualize vessel labels
cmd = {src.fs};
cmd{end+1} = 'freeview \';
cmd{end+1} = ['-v ' mriSeg.fspec ':resample=nearest \'];
cmd{end+1} = ['-v ' mriTof.fspec ':resample=nearest \'];
cmd{end+1} = ['-v ' vfMRI ':resample=nearest \'];
cmd{end+1} = ['-v ' mriLabels.fspec ':resample=nearest'];
disp(strjoin(cmd,newline)); clear cmd;
% select vessels
vesselIdx = sort([34 18 7 8 25 22 6 17 14 16 80 19 10 13 33 21 11 15]); %1->supuriously connected to the sagittal sinus
% rewrite labels with only selected vessels
mriLabels.fspec = replace(mriLabels.fspec, 'tofLabels1.nii', 'tofLabels2.nii');
if forceThis || ~exist(mriLabels.fspec,'file')
    mriLabels.vol(~ismember(mriLabels.vol,vesselIdx)) = 0;
    MRIwrite(mriLabels, mriLabels.fspec);
else
    mriLabels = MRIread(mriLabels.fspec);
end
% split into individual vessels
clear vessels
for v = 1:length(vesselIdx)
    cmd = {src.fs};
    vessels(v).label      = vesselIdx(v);
    vessels(v).ref.fLabel = mriLabels.fspec;
    vessels(v).ref.fTof   = mriTof.fspec;
    % vessels(v).fTof     = mriTof.fspec;
    % vessels(v).fTofList = tofList;
    % vessels(v).fSegList = tofSegList;

    % find bounding box for orginal resolution
    mri = mriLabels;
    mri.vol = mri.vol==vesselIdx(v);
    [dim1, dim2, dim3] = ind2sub(size(mri.vol), find(mri.vol));
    dim1 = [min(dim1) max(dim1)]; dim1 = dim1 + [-1 1]; dim1(dim1<1) = 1;
    dim2 = [min(dim2) max(dim2)]; dim2 = dim2 + [-1 1]; dim2(dim2<1) = 1;
    dim3 = [min(dim3) max(dim3)]; dim3 = dim3 + [-1 1]; dim3(dim3<1) = 1;
    dim1(dim1>size(mri.vol,1)) = size(mri.vol,1);
    dim2(dim2>size(mri.vol,2)) = size(mri.vol,2);
    dim3(dim3>size(mri.vol,3)) = size(mri.vol,3);
    if mod(mean(dim1),1); if dim1(1) ~= 1; dim1(1) = dim1(1) - 1; else; dim1(2) = dim1(2) + 1; end; end
    if mod(mean(dim2),1); if dim2(1) ~= 1; dim2(1) = dim2(1) - 1; else; dim2(2) = dim2(2) + 1; end; end
    if mod(mean(dim3),1); if dim3(1) ~= 1; dim3(1) = dim3(1) - 1; else; dim3(2) = dim3(2) + 1; end; end
    % crop original resolution label to a reference grid
    vessels(v).cropRef.center = [ mean(dim2)  mean(dim1)          1 ]-1;
    vessels(v).cropRef.size   = [range(dim2) range(dim1) range(dim3)]+1;
    vessels(v).cropRef.fLabel = replace(vessels(v).ref.fLabel, 'tofLabels2.nii', ['refLabel-vessel' sprintf('%02d',vesselIdx(v)) '.nii.gz']);
    if forceThis || ~exist(vessels(v).cropRef.fLabel,'file')
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cmd{end+1} = ['--crop '       strjoin(arrayfun(@num2str, vessels(v).cropRef.center, 'UniformOutput', false), ' ') ' \'];
        cmd{end+1} = ['--cropsize '   strjoin(arrayfun(@num2str, vessels(v).cropRef.size  , 'UniformOutput', false), ' ') ' \'];
        cmd{end+1} = [vessels(v).ref.fLabel ' \'];
        cmd{end+1} = vessels(v).cropRef.fLabel;
        system(strjoin(cmd,newline),'-echo');
    end
    % turn it into an all-one mask for easier debuging
    vessels(v).cropRef.fMask  = replace(vessels(v).cropRef.fLabel,'/refLabel-','/refMask-');
    if forceThis || ~exist(vessels(v).cropRef.fMask,'file')
        tmp = MRIread(vessels(v).cropRef.fLabel);
        tmp.vol(:) = uint8(1);
        MRIwrite(tmp, vessels(v).cropRef.fMask);
    end
    

    
    % resample reference grids (label and mask) to all resolutions
    vessels(v).cropRef.usList = usList;
    cmd = {src.fs};
    for i = 1:length(vessels(v).cropRef.usList)
        vessels(v).cropRef.fLabelList{i} = fullfile(fileparts(vessels(v).cropRef.fLabel), ['vessel-' sprintf('%02d',vesselIdx(v)) '_upsampledX-' num2str(usList(i)) '_label.nii.gz']);
        vessels(v).cropRef.fMaskList{i}  = fullfile(fileparts(vessels(v).cropRef.fLabel), ['vessel-' sprintf('%02d',vesselIdx(v)) '_upsampledX-' num2str(usList(i)) '_mask.nii.gz']);
        if forceThis || ~exist(vessels(v).cropRef.fLabelList{i},'file')
            cmd{end+1} = 'mri_convert \';
            cmd{end+1} = ['--resample_type nearest \'];
            cmd{end+1} = ['--upsample ' num2str(usList(i)) ' \'];
            cmd{end+1} = [vessels(v).cropRef.fLabel ' \'];
            cmd{end+1} = vessels(v).cropRef.fLabelList{i};
        end
        if forceThis || ~exist(vessels(v).cropRef.fMaskList{i},'file')
            cmd{end+1} = 'mri_convert \';
            cmd{end+1} = ['--resample_type nearest \'];
            cmd{end+1} = ['--upsample ' num2str(usList(i)) ' \'];
            cmd{end+1} = [vessels(v).cropRef.fMask ' \'];
            cmd{end+1} = vessels(v).cropRef.fMaskList{i};
        end
    end
    if length(cmd) > 1
        system(strjoin(cmd,newline),'-echo');
    end
    clear cmd;


    cmd = {src.fs};
    % crop (resample --like) every vesselboost outputs
    for p = 1:length(tofSegList)
        % find which upsampled version
        usIdx = strsplit(tofSegList{p},'Ups'); usIdx = strsplit(usIdx{end},'x'); usIdx = str2num(usIdx{1}); usIdx = find(vessels(v).cropRef.usList == usIdx);
        vessels(v).upsampleFactor(p) = vessels(v).cropRef.usList(usIdx);
        % seg
        in = tofSegList{p};
        vessels(v).fSeg{p} = fullfile(fileparts(in), ['vessel-' sprintf('%02d',vesselIdx(v)) '_seg.nii.gz']);
        if forceThis || ~exist(vessels(v).fSeg{p},'file')
            cmd{end+1} = 'mri_convert \';
            cmd{end+1} = ['--resample_type nearest \'];
            cmd{end+1} = ['--like ' vessels(v).cropRef.fMaskList{usIdx} ' \'];
            cmd{end+1} = [in ' \'];
            cmd{end+1} = [vessels(v).fSeg{p}];
        end
        % tof
        in = tofList{p};
        vessels(v).fTof{p} = fullfile(fileparts(in), ['vessel-' sprintf('%02d',vesselIdx(v)) '_tof.nii.gz']);
        if forceThis || ~exist(vessels(v).fTof{p},'file')
            cmd{end+1} = 'mri_convert \';
            cmd{end+1} = ['--resample_type nearest \'];
            cmd{end+1} = ['--like ' vessels(v).cropRef.fMaskList{usIdx} ' \'];
            cmd{end+1} = [in ' \'];
            cmd{end+1} = [vessels(v).fTof{p}];
        end
    end


    % run commands
    disp('----------------------------------');
    disp('----------------------------------');
    disp(['Vessel ' sprintf('%02d',vesselIdx(v)) ' (' num2str(v) '/' num2str(length(vesselIdx)) ') writing...']);
    disp('----------------------------------');
    disp('----------------------------------');
    if length(cmd) > 1
        tic;
        system(strjoin(cmd,newline),'-echo');
        disp('----------------------------------');
        disp('----------------------------------');
        disp(['Vessel ' sprintf('%02d',vesselIdx(v)) ' (' num2str(v) '/' num2str(length(vesselIdx)) ') writing done...']);
        toc;
        disp('----------------------------------');
        disp('----------------------------------');
    else
        disp('----------------------------------');
        disp('----------------------------------');
        disp(['Vessel ' sprintf('%02d',vesselIdx(v)) ' (' num2str(v) '/' num2str(length(vesselIdx)) ') writing alreadydone...']);
        disp('----------------------------------');
        disp('----------------------------------');
    end
    clear cmd;

end
%% %%%%%%%%%%%%%%%%%%%%%
info.vessel.pointer = vessels;
info.vessel.pointerFile = fullfile(info.project.code,mfilename,['sub-' num2str(S) '_vesselPointer' datestr(now,'yyyymmddHHMMSS') '.mat']);
if ~exist(fileparts(info.vessel.pointerFile),'dir'); mkdir(fileparts(info.vessel.pointerFile)); end
if rxivThis || ~exist(info.vessel.pointerFile,'file')
    save(info.vessel.pointerFile,'info');
    disp(['Vessel pointer saved to:' newline info.vessel.pointerFile]);
end


