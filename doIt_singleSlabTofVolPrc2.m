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


forceThis = 1;
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofCrop;
tofCropPrc;
tofCropPrcSeg;




forceThis = 1;
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
    % Upsample (cubic interpolation) preprocessed tof
    tofCropPrcScl{i} = replace(tofCropPrc, '/cropPrc/', ['/cropPrcScale' num2str(usList(i)) '/']);
    if ~exist(fileparts(tofCropPrcScl{i}),'dir'); mkdir(fileparts(tofCropPrcScl{i})); end
    if forceThis || ~exist(tofCropPrcScl{i},'file')
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cmd{end+1} = ['--voxsize ' strjoin(arrayfun(@(x) num2str(x,'%.15g'), res./usList(i), 'UniformOutput', false), ' ') ' \'];
        cmd{end+1} = ['--resample_type cubic \'];
        cmd{end+1} = [tofCropPrc ' \'];
        cmd{end+1} = [tofCropPrcScl{i}];
        system(strjoin(cmd,newline),'-echo');
    end

    % Vesselboost
    tofCropPrcSclSeg{i} = replace(tofCropPrcScl{i},'/tof.nii.gz','/seg.nii.gz');
    vesselboost_prediction(...
        tofCropPrcScl{i}   ,...
        tofCropPrcSclSeg{i},...
        [],...
        vesselBoostModel,4);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofCropPrcScl;
tofCropPrcSclSeg;





forceThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%
%% Individualize vessels
%%%%%%%%%%%%%%%%%%%%%%%%
% find connected components (vessels) based on segmentation from non-upsampled but preprocessed tof
mriSeg     = MRIread(tofCropPrcSeg);
mriSeg.vol = logical(mriSeg.vol);
mriTof     = MRIread(tofCropPrc);
mriLabels  = mriSeg;
mriLabels.fspec = replace(mriSeg.fspec,'/seg.nii.gz','/label.nii.gz');
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
mriLabels.fspec = replace(mriLabels.fspec, 'label.nii.gz', 'selectedLabel.nii.gz');
if forceThis || ~exist(mriLabels.fspec,'file')
    mriLabels.vol(~ismember(mriLabels.vol,vesselIdx)) = 0;
    MRIwrite(mriLabels, mriLabels.fspec);
else
    mriLabels = MRIread(mriLabels.fspec);
end


% Crop out individual vessels
clear vessels
for v = 1:length(vesselIdx)
    disp(' ')
    disp('----------------------------------');
    disp('----------------------------------');
    disp(['Vessel ' sprintf('%02d',vesselIdx(v)) ' (' num2str(v) '/' num2str(length(vesselIdx)) ') cropping out...']);
    disp('----------------------------------');
    disp('----------------------------------');

    cmd = {src.fs};
    vessels(v).label      = vesselIdx(v);
    vessels(v).ref.fLabel = mriLabels.fspec;
    vessels(v).ref.fTof   = mriTof.fspec;
    
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
    vessels(v).cropRef.fMask  = replace(vessels(v).ref.fLabel, 'selectedLabel.nii.gz', ['vessel-' sprintf('%02d',vesselIdx(v)) '_mask.nii.gz']);
    if forceThis || ~exist(vessels(v).cropRef.fMask,'file')
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cmd{end+1} = ['--crop '       strjoin(arrayfun(@num2str, vessels(v).cropRef.center, 'UniformOutput', false), ' ') ' \'];
        cmd{end+1} = ['--cropsize '   strjoin(arrayfun(@num2str, vessels(v).cropRef.size  , 'UniformOutput', false), ' ') ' \'];
        cmd{end+1} = [vessels(v).ref.fLabel ' \'];
        cmd{end+1} = vessels(v).cropRef.fMask;
        system(strjoin(cmd,newline),'-echo');
        mri = MRIread(vessels(v).cropRef.fMask);
        mri.vol = mri.vol==vessels(v).label;
        MRIwrite(mri, vessels(v).cropRef.fMask);
    end

    vessels(v).cropRef.fMask;
    vessels(v).cropRef.res       = res;
    vessels(v).cropRef.usList    = usList;
    vessels(v).cropRef.resList   = cell(1, length(vessels(v).cropRef.usList));
    vessels(v).cropRef.fMaskList = cell(1, length(vessels(v).cropRef.usList));
    vessels(v).fMaskList = cell(length(vessels(v).cropRef.usList), 1);
    for p = 1:length(vessels(v).cropRef.usList)
        disp(' ')
        disp(['---------------------------']);
        disp(['scale=' num2str(vessels(v).cropRef.usList(p)) ' (' num2str(p) '/' num2str(length(vessels(v).cropRef.usList)) ') processing...']);
        disp(['---------------------------']);

        % To obtain a reference grid (--like), upsample the vessel mask derived from original resolution data
        vessels(v).cropRef.resList{p}   = res./vessels(v).cropRef.usList(p);
        vessels(v).cropRef.fMaskList{p} = replace(vessels(v).cropRef.fMask, 'mask.nii.gz', ['maskScale' num2str(usList(p)) '.nii.gz']);
        if forceThis || ~exist(vessels(v).cropRef.fMaskList{p},'file')
            cmd = {src.fs};
            cmd{end+1} = 'mri_convert \';
            cmd{end+1} = ['--resample_type nearest \'];
            cmd{end+1} = ['--voxsize ' strjoin(arrayfun(@(x) num2str(x,'%.15g'), vessels(v).cropRef.resList{p}, 'UniformOutput', false), ' ') ' \'];
            cmd{end+1} = [vessels(v).cropRef.fMask ' \'];
            cmd{end+1} =  vessels(v).cropRef.fMaskList{p};
            system(strjoin(cmd,newline),'-echo');
        end
        % Use the reference grid to crop segmentation from all scales, then upsample (nearest neighbor) to a mask at the highest scale
        vessels(v).fMaskList{p} = fullfile(fileparts(tofCropPrcSclSeg{p}), ['vessel-' sprintf('%02d',vesselIdx(v)) '_mask.nii.gz']);
        if forceThis || ~exist(vessels(v).fMaskList{p},'file')
            % crop (resample --like) vesselboost segmentation
            cmd = {src.fs};
            cmd{end+1} = 'mri_convert \';
            cmd{end+1} = ['--resample_type nearest \'];
            cmd{end+1} = ['--like ' vessels(v).cropRef.fMaskList{p} ' \'];
            cmd{end+1} = [tofCropPrcSclSeg{p} ' \'];
            cmd{end+1} = vessels(v).fMaskList{p};
            system(strjoin(cmd,newline),'-echo');
            % match candidate vessels in seg from upsampled data to the one vessel in the mask derived from original resolution data
            mriRef = MRIread(vessels(v).cropRef.fMaskList{p});
            mri    = MRIread(vessels(v).fMaskList{p});
            CCRef = bwconncomp(mriRef.vol,26);
            CC    = bwconncomp(mri.vol   ,26);
            n = zeros(length(CC.PixelIdxList),1);
            for c = 1:length(CC.PixelIdxList)
                n(c) = nnz(ismember(CCRef.PixelIdxList{1},CC.PixelIdxList{c}));
            end
            [~,b] = max(n);
            % turn the segmentation into a mask of the best-match vessel
            eMask = true(size(mri.vol));
            eMask(CC.PixelIdxList{b}) = false;
            mri.vol(eMask) = false;
            MRIwrite(mri, vessels(v).fMaskList{p});
            % upsample to the highest scale
            [~,b] = max(vessels(v).cropRef.usList);
            maxRes = vessels(v).cropRef.resList{b}
            if any(vessels(v).cropRef.resList{p}~=maxRes)
                cmd = {src.fs};
                cmd{end+1} = 'mri_convert \';
                cmd{end+1} = ['--resample_type nearest \'];
                cmd{end+1} = ['--voxsize ' strjoin(arrayfun(@(x) num2str(x,'%.15g'), maxRes, 'UniformOutput', false), ' ') ' \'];
                cmd{end+1} = [vessels(v).fMaskList{p} ' \'];
                cmd{end+1} = vessels(v).fMaskList{p};
                system(strjoin(cmd,newline),'-echo');
            end
        end
        disp(['---------------------------']);
        disp(['scale=' num2str(vessels(v).cropRef.usList(p)) ' (' num2str(p) '/' num2str(length(vessels(v).cropRef.usList)) ') done.']);
        disp(['---------------------------']);
        disp(' ')
    end
    disp('----------------------------------');
    disp('----------------------------------');
    disp(['Vessel ' sprintf('%02d',vesselIdx(v)) ' (' num2str(v) '/' num2str(length(vesselIdx)) ') cropping out done.']);
    disp('----------------------------------');
    disp('----------------------------------');
    disp(' ')
end
%% %%%%%%%%%%%%%%%%%%%%%
vessels;


forceThis = 1;
%%%%%%%%%%%%%%%%%%%%%%
%% Save vessel pointer
%%%%%%%%%%%%%%%%%%%%%%
info.vessel.pointer = vessels;
info.vessel.pointerFile = fullfile(info.project.code,mfilename,['sub-' num2str(S) '_vesselPointer' datestr(now,'yyyymmddHHMMSS') '.mat']);
if ~exist(fileparts(info.vessel.pointerFile),'dir'); mkdir(fileparts(info.vessel.pointerFile)); end
if forceThis || ~exist(info.vessel.pointerFile,'file')
    save(info.vessel.pointerFile,'info');
    disp(['Vessel pointer saved to:' newline info.vessel.pointerFile]);
end


