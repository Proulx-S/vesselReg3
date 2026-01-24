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

S=2; arxvDir = '/local/users/Proulx-S/db/vsmDiamCenSur/sub-vsmDiamCenSurP2_acq-vfMRI_prsc-dflt_venc-none';
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
            in = fullfile(arxvDir,'tofBrain.nii.gz');
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
arxvThis  = 0;
usList = [1 2 3 4 6 8];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vesselboost on cropped -> upsampled -> bias-field-corrected -> denoised tof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    arxvFile = strsplit(tofSegList1{i}, '/'); arxvFile = fullfile(arxvDir,[arxvFile{end-1} '.nii.gz']);
    if arxvThis
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cmd{end+1} = [tofSegList1{i} ' \'];
        cmd{end+1} = arxvFile;
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
    arxvFile = strsplit(tofSegList2{i}, '/'); arxvFile = fullfile(arxvDir,[arxvFile{end-1} '.nii.gz']);
    if arxvThis
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cmd{end+1} = [tofSegList2{i} ' \'];
        cmd{end+1} = arxvFile;
        system(strjoin(cmd,newline),'-echo');
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(tofSegList1)
disp(tofSegList2)


forceThis = 0;
arxvThis  = 0;
usList = [1 2 3 4 6 8];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vesselboost on cropped -> bias-field-corrected -> denoised -> upsampled tof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    arxvFile = strsplit(tofSegList3{i}, '/'); arxvFile = fullfile(arxvDir,[arxvFile{end-1} '.nii.gz']);
    if arxvThis
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cmd{end+1} = [tofSegList3{i} ' \'];
        cmd{end+1} = arxvFile;
        system(strjoin(cmd,newline),'-echo');
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(tofSegList3)



forceThis = 1;
arxvThis  = 1;
%%%%%%%%%%%%%%%%%%%%%%%%
%% Individualize vessels
%%%%%%%%%%%%%%%%%%%%%%%%
% find connected components (vessels)
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
disp(strjoin(cmd,newline));
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
% write mask for each selected vessel
clear vessels
for v = 1:length(vesselIdx)
    vessels(v).fname = replace(mriLabels.fspec, 'tofLabels2.nii', ['tofVessel' sprintf('%02d',vesselIdx(v)) '.nii']);
    % create mask
    mri = MRIread(mriLabels.fspec);
    mri.vol = mri.vol==vesselIdx(v);
    mri.fspec = vessels(v).fname;
    % find bounding box
    [dim1, dim2, dim3] = ind2sub(size(mri.vol), find(mri.vol));
    dim1 = [min(dim1) max(dim1)]; dim1 = dim1 + [-1 1]; dim1(dim1<1) = 1; dim1(dim1>size(mriLabels.vol,1)) = size(mriLabels.vol,1);
    dim2 = [min(dim2) max(dim2)]; dim2 = dim2 + [-1 1]; dim2(dim2<1) = 1; dim2(dim2>size(mriLabels.vol,2)) = size(mriLabels.vol,2);
    dim3 = [min(dim3) max(dim3)]; dim3 = dim3 + [-1 1]; dim3(dim3<1) = 1; dim3(dim3>size(mriLabels.vol,3)) = size(mriLabels.vol,3);
    if mod(mean(dim1),1); if dim1(1) ~= 1; dim1(1) = dim1(1) - 1; else; dim1(2) = dim1(2) + 1; end; end
    if mod(mean(dim2),1); if dim2(1) ~= 1; dim2(1) = dim2(1) - 1; else; dim2(2) = dim2(2) + 1; end; end
    if mod(mean(dim3),1); if dim3(1) ~= 1; dim3(1) = dim3(1) - 1; else; dim3(2) = dim3(2) + 1; end; end
    vessels(v).dim1 = dim1;
    vessels(v).dim2 = dim2;
    vessels(v).dim3 = dim3;
    % crop mask
    if forceThis || ~exist(vessels(v).fname,'file')
        MRIwrite(mri, mri.fspec);
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cropCenter = [mean(dim2) mean(dim1) 0];
        cropSize = [range(dim2)+1 range(dim1)+1 range(dim3)+1];
        cmd{end+1} = ['--crop '     strjoin(arrayfun(@num2str, cropCenter-1, 'UniformOutput', false), ' ') ' \'];
        cmd{end+1} = ['--cropsize ' strjoin(arrayfun(@num2str, cropSize  , 'UniformOutput', false), ' ') ' \'];
        cmd{end+1} = ['--slice-crop ' strjoin(arrayfun(@num2str, dim3-1, 'UniformOutput', false), ' ') ' \'];
        cmd{end+1} = [vessels(v).fname ' \'];
        cmd{end+1} = [vessels(v).fname];
        system(strjoin(cmd,newline),'-echo');
    end
    % confirm proper cropping
    if 0
        mri1 = MRIread(vessels(v).fname);
        mri1.vol = mri.vol(dim1(1):dim1(2), dim2(1):dim2(2), dim3(1):dim3(2));
        mri1.fspec = replace(vessels(v).fname, '.nii', 'X.nii');
        MRIwrite(mri1, mri1.fspec);
        mri2 = MRIread(vessels(v).fname);
        size(mri1.vol)
        size(mri2.vol)
        max(abs(mri1.vol(:) - mri2.vol(:)))/max(abs([mri1.vol(:); mri2.vol(:)]))*100
        cmd = {src.fs};
        cmd{end+1} = 'freeview \';
        cmd{end+1} = [mriTof.fspec ' \'];
        cmd{end+1} = [mri1.fspec ' \'];
        cmd{end+1} = [mri2.fspec];
        disp(strjoin(cmd,newline));
    end
    % archive
    arxvFile = strsplit(vessels(v).fname, '/'); arxvFile = fullfile(arxvDir,arxvFile{end}); arxvFile = replace(arxvFile, '.nii', '.nii.gz');
    if arxvThis
        cmd = {src.fs};
        cmd{end+1} = 'mri_convert \';
        cmd{end+1} = [vessels(v).fname ' \'];
        cmd{end+1} = arxvFile;
        system(strjoin(cmd,newline),'-echo');
    end
end
%% %%%%%%%%%%%%%%%%%%%%%
disp(vessels)


return















I1 = 4;
tof_upsampled = replace(tof, '/tof/', ['/tofUps' num2str(I1) 'x/']);
if ~exist(fileparts(tof_upsampled),'dir'); mkdir(fileparts(tof_upsampled)); end
if forceThis || ~exist(tof_upsampled,'file')
    cmd = {src.fs};
    cmd{end+1} = 'mri_convert \';
    cmd{end+1} = ['--resample_type cubic \'];
    cmd{end+1} = ['--upsample ' num2str(I1) ' \'];
    cmd{end+1} = [tof ' \'];
    cmd{end+1} = [tof_upsampled];
    system(strjoin(cmd,newline),'-echo');
end
tofOrig = tof;
tof = tof_upsampled;





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tofSegList1
tofSegList2
tofSegList3




return






% Skullstrip (extract the brain) using FreeSurfer's machine learning approach: mri_synthstrip
% Output: 'tof_brain.nii.gz' in the same folder as 'tof'
tof_brain     = replace(tof,'/tof/','/tofBrain/'); if ~exist(fileparts(tof_brain),'dir'); mkdir(fileparts(tof_brain)); end
tof_brainMask = replace(tof,'/tof/','/tofBrainMask/'); if ~exist(fileparts(tof_brainMask),'dir'); mkdir(fileparts(tof_brainMask)); end
tof_brainDist = replace(tof,'/tof/','/tofBrainDist/'); if ~exist(fileparts(tof_brainDist),'dir'); mkdir(fileparts(tof_brainDist)); end

cmd = {src.fs};
cmd{end+1} = 'mri_synthstrip \';
cmd{end+1} = ['-i ' tof ' \'];
cmd{end+1} = ['-m ' tof_brainMask ' \'];
cmd{end+1} = ['-d ' tof_brainDist ' \'];
cmd{end+1} = ['-o ' tof_brain];
disp(strjoin(cmd,newline))
system(strjoin(cmd,newline),'-echo');



return
% tof = tofOrig;
forceThis = 0;
I1 = 4;
tof_upsampled = replace(tof, '/tof/', ['/tofUps' num2str(I1) 'x/']);
if ~exist(fileparts(tof_upsampled),'dir'); mkdir(fileparts(tof_upsampled)); end
if forceThis || ~exist(tof_upsampled,'file')
    cmd = {src.fs};
    cmd{end+1} = 'mri_convert \';
    cmd{end+1} = ['--resample_type cubic \'];
    cmd{end+1} = ['--upsample ' num2str(I1) ' \'];
    cmd{end+1} = [tof ' \'];
    cmd{end+1} = [tof_upsampled];
    system(strjoin(cmd,newline),'-echo');
end
tofOrig = tof;
tof = tof_upsampled;




forceThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vesselboost vessel segmentation of tof volume data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tof_vesselSeg = fullfile(projectScratch,'tof','seg'); if ~exist(tof_vesselSeg,'dir'); mkdir(tof_vesselSeg); end
tof_vesselSeg = fullfile(tof_vesselSeg,'tof.nii.gz');
ps_path = fullfile(projectScratch,'ps'); if ~exist(ps_path,'dir'); mkdir(ps_path); end
if forceThis || ~exist(tof_vesselSeg,'file')
    vesselboost_prediction(fileparts(tof),fileparts(tof_vesselSeg),ps_path,vesselBoostModel,3); % should also try with prepmode=3 (bias field correction and denoising preprocessing)
end
%reduce precision and don't compress
tmp = replace(tof_vesselSeg, '.nii.gz', '.nii');
cmd = {src.fs};
cmd{end+1} = 'mri_convert \';
% cmd{end+1} = '--out_data_type int \';
cmd{end+1} = '--out_data_type uchar \';
cmd{end+1} = [tof_vesselSeg ' \'];
cmd{end+1} = tmp;
system(strjoin(cmd,newline),'-echo');
tof_vesselSeg = tmp;


% %upsample
% I = 8;
% cmd = {src.fs};
% cmd{end+1} = 'mri_convert \';
% cmd{end+1} = ['--upsample ' num2str(I) ' \'];
% cmd{end+1} = ['--resample_type nearest \'];
% cmd{end+1} = [tof_vesselSegX ' \'];
% cmd{end+1} = tof_vesselSegX;
% system(strjoin(cmd,newline),'-echo');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Individualize and skeletonize vessels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find connected components (vessels)
mriSeg = MRIread(tof_vesselSeg);
mriSeg.vol = logical(mriSeg.vol);
disp('bwconncomp: running...');
CC = bwconncomp(mriSeg.vol,26);
disp('bwconncomp: done');
% sort vessels from big to small vessel network volume
[~,b] = sort(cellfun('length',CC.PixelIdxList),'descend');
CC.PixelIdxList = CC.PixelIdxList(b);
% label indivudal vessels in a single map for visualization
mriVesselLabel = mriSeg;
mriVesselLabel.vol = zeros(size(mriSeg.vol));
for v = 1:length(CC.PixelIdxList)
    mriVesselLabel.vol(CC.PixelIdxList{v}) = v;
end
mriVesselLabel.fspec = replace(tof_vesselSeg,'.nii','_vesselLabels.nii');
MRIwrite(mriVesselLabel, mriVesselLabel.fspec);
% visualize vessel labels
cmd = {src.fs};
cmd{end+1} = 'freeview \';
cmd{end+1} = ['-v ' tof ' \'];
cmd{end+1} = ['-v ' vfMRI ' \'];
cmd{end+1} = ['-v ' tof_vesselSeg ' \'];
cmd{end+1} = ['-v ' mriVesselLabel.fspec];
disp(strjoin(cmd,newline));
% select vessels
% vesselIdx = sort([33 17 19 20 13 52 15 10]); % 1 is supuriously connected to the sagittal sinus
vesselIdx = sort([31 19 18 123 10 17 8 11]); % vesselboos on upsampled tof x 2 and prepmode=3 % 2 and 3 are supuriously connected to sinus
return


% process selected vessels --method 1
mriTof = MRIread(tof);
I = 8/I1; P = I/2;
for v = 1:length(vesselIdx)
    %crop
    mriVessel = mriVesselLabel;
    mriVessel.fspec = replace(mriVesselLabel.fspec, '.nii', '_tmpVesselOrig.nii');
    mriVessel.vol = mriVesselLabel.vol==vesselIdx(v);
    inds = find(mriVessel.vol);
    [dim1, dim2, dim3] = ind2sub(size(mriVessel.vol), inds);
    mriVessel.bbox.min = [min(dim1), min(dim2), min(dim3)];
    mriVessel.bbox.max = [max(dim1), max(dim2), max(dim3)];
    mriVessel.bbox.size = mriVessel.bbox.max - mriVessel.bbox.min + 1;
    dim1S = mriVessel.bbox.min(1); dim1E = mriVessel.bbox.max(1); dim2S = mriVessel.bbox.min(2); dim2E = mriVessel.bbox.max(2); dim3S = mriVessel.bbox.min(3); dim3E = mriVessel.bbox.max(3);
    mriVessel.vol = mriVessel.vol(dim1S:dim1E, dim2S:dim2E, dim3S:dim3E);
    mriVessel.volTof = mriTof.vol(dim1S:dim1E, dim2S:dim2E, dim3S:dim3E);
    %zero pad
    mriVessel.vol = padarray(mriVessel.vol,[P P P],0,'both');
    mriVessel.volTof = padarray(mriVessel.volTof,[P P P],0,'both');
    %upsample
    mriVessel.fspec = replace(mriVesselLabel.fspec, '.nii', '_tmpVesselProcessed.nii');
    mriVessel.vol = imresize3(double(mriVessel.vol),I,'lanczos3')>0.5;
    mriVessel.volTof = imresize3(double(mriVessel.volTof),I,'lanczos3');
    %closing (dilate then erode) fills gaps
    mriVessel.vol = imclose(mriVessel.vol,strel('sphere',I));
    %opening (erode then dilate) removes corner protrusions and small artifacts
    mriVessel.vol = imopen(mriVessel.vol,strel('sphere',1));
    %smooth
    mriVessel.vol = imgaussfilt3(double(mriVessel.vol),I)>1/I;
    %skeletonize
    mriVessel.vol = bwskel(mriVessel.vol,'MinBranchLength',I*2);
    %threshold distance transform
    mriVessel.vol = bwdist(mriVessel.vol)<I;
    MRIwrite(mriVessel, mriVessel.fspec); disp(mriVessel.fspec)
    %surface
    surfFile = replace(mriVessel.fspec, '.nii', '.vtk');
    vmtk_surfFromSeg(mriVessel.fspec, surfFile)
    %centerline
    lineFile = replace(surfFile, '.vtk', '.line.vtk');
    cmd = vmtk_centerlinesFromSurf(surfFile, lineFile,[],[],1);
    disp(strjoin(cmd,newline));

    %write vessel tof
    mri = mriVessel; mri.vol = mri.volTof; mri.fspec = replace(mriVessel.fspec, '.nii', '_tmpVesselTof.nii'); MRIwrite(mri, mri.fspec); disp(mri.fspec)
    %view
    vmtk_viewVolAndSurf(mri.fspec, lineFile)
end


% process selected vessels --method 2
mriTof = MRIread(tof);
I = 8/I1; P = I/2;
for v = 1:length(vesselIdx)
    %crop
    mriVessel = mriVesselLabel;
    mriVessel.fspec = replace(mriVesselLabel.fspec, '.nii', '_tmpVesselOrig.nii');
    mriVessel.vol = mriVesselLabel.vol==vesselIdx(v);
    inds = find(mriVessel.vol);
    [dim1, dim2, dim3] = ind2sub(size(mriVessel.vol), inds);
    mriVessel.bbox.min = [min(dim1), min(dim2), min(dim3)];
    mriVessel.bbox.max = [max(dim1), max(dim2), max(dim3)];
    mriVessel.bbox.size = mriVessel.bbox.max - mriVessel.bbox.min + 1;
    dim1S = mriVessel.bbox.min(1); dim1E = mriVessel.bbox.max(1); dim2S = mriVessel.bbox.min(2); dim2E = mriVessel.bbox.max(2); dim3S = mriVessel.bbox.min(3); dim3E = mriVessel.bbox.max(3);
    mriVessel.vol = mriVessel.vol(dim1S:dim1E, dim2S:dim2E, dim3S:dim3E);
    mriVessel.volTof = mriTof.vol(dim1S:dim1E, dim2S:dim2E, dim3S:dim3E);
    %zero pad
    mriVessel.vol = padarray(mriVessel.vol,[P P P],0,'both');
    mriVessel.volTof = padarray(mriVessel.volTof,[P P P],0,'both');
    %upsample
    mriVessel.fspec = replace(mriVesselLabel.fspec, '.nii', '_tmpVesselProcessed.nii');
    mriVessel.vol = imresize3(double(mriVessel.vol),I,'lanczos3')>0.75;
    mriVessel.volTof = imresize3(double(mriVessel.volTof),I,'lanczos3');
    %skeletonize
    mriVessel.vol = bwskel(mriVessel.vol,'MinBranchLength',I*2);
    %threshold distance transform
    mriVessel.vol = bwdist(mriVessel.vol)<2;
    MRIwrite(mriVessel, mriVessel.fspec); disp(mriVessel.fspec)



    %closing (dilate then erode) fills gaps
    mriVessel.vol = imclose(mriVessel.vol,strel('sphere',I));
    %opening (erode then dilate) removes corner protrusions and small artifacts
    mriVessel.vol = imopen(mriVessel.vol,strel('sphere',1));



    %smooth
    mriVessel.vol = imgaussfilt3(double(mriVessel.vol),I)>1/I;
    %surface
    surfFile = replace(mriVessel.fspec, '.nii', '.vtk');
    vmtk_surfFromSeg(mriVessel.fspec, surfFile)
    %centerline
    lineFile = replace(surfFile, '.vtk', '.line.vtk');
    cmd = vmtk_centerlinesFromSurf(surfFile, lineFile,[],[],1);
    disp(strjoin(cmd,newline));

    %write vessel tof
    mri = mriVessel; mri.vol = mri.volTof; mri.fspec = replace(mriVessel.fspec, '.nii', '_tmpVesselTof.nii'); MRIwrite(mri, mri.fspec); disp(mri.fspec)
    %view
    vmtk_viewVolAndSurf(mri.fspec, lineFile)
end



output_vtk = skeleton_to_graph_vtk(skeleton_nii, output_vtk, connectivity, forceThis)


return


disp(['bwskel: done']);
% sort vessels from long to short vessel network
[~,b] = sort(nx,'descend');
mrixVesselSeg  = mrixVesselSeg(b);
mrixVesselSkel = mrixVesselSkel(b);
nx             = nx(b);
% label indivudal vessels in a single (downsampled) map for visualization
mriSeg = MRIread(tof_vesselSeg);
mriVesselLabel = mriSeg;
mriVesselLabel.vol = zeros(size(mriSeg.vol), 'uint8');
for v = 1:length(mrixVesselSeg)
    mriVesselLabel.vol(imresize3(mrixVesselSeg{v}.vol,1/I,'nearest')) = uint8(v);
    mrixVesselSeg{v}.fspec  = replace(mrixVesselSeg{v}.fspec, '.nii', ['_vessel' sprintf('%03d',v) '.nii']);
    mrixVesselSkel{v}.fspec = replace(mrixVesselSkel{v}.fspec, '.nii', ['_vessel' sprintf('%03d',v) '_skel.nii']);
end
mriVesselLabel.fspec = replace(tof_vesselSeg,'.nii.gz','_vesselLabels.nii.gz');
MRIwrite(mriVesselLabel, mriVesselLabel.fspec);
% visualize and select vessels
cmd = {src.fs};
cmd{end+1} = 'freeview \';
cmd{end+1} = ['-v ' tof ' \'];
cmd{end+1} = ['-v ' vfMRI ' \'];
cmd{end+1} = ['-v ' tof_vesselSeg ' \'];
cmd{end+1} = ['-v ' mriVesselLabel.fspec];
disp(strjoin(cmd,newline));
vesselIdx = [6 35 20 18 19 51 15 13 21 32 8 9 29 11 39 108 24 50]; % 2 is supuriously connected to the sagittal sinus
mrixVesselSeg = mrixVesselSeg(vesselIdx);
mrixVesselSkel = mrixVesselSkel(vesselIdx);
% % write selected individual vessel segmentation and skeleton masks
% for v = 1:length(vesselIdx)
%     MRIwrite(mriVesselSeg{vesselIdx(v)} ,mriVesselSeg{vesselIdx(v)}.fspec);
%     MRIwrite(mriVesselSkel{vesselIdx(v)},mriVesselSkel{vesselIdx(v)}.fspec);
% end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Crop individual vessel masks and dilate/erode to fill gaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vesselSurf = cell(length(mrixVesselSeg), 1);
for v = 1:length(mrixVesselSeg)
    %crop
    vesselCrop = mrixVesselSeg{v}.vol;
    inds = find(vesselCrop);
    [y, x, z] = ind2sub(size(vesselCrop), inds);
    bbox.min = [min(x), min(y), min(z)];
    bbox.max = [max(x), max(y), max(z)];
    bbox.size = bbox.max - bbox.min + 1;
    xS = bbox.min(1); xE = bbox.max(1); yS = bbox.min(2); yE = bbox.max(2); zS = bbox.min(3); zE = bbox.max(3);
    vesselCrop = vesselCrop(xS:xE, yS:yE, zS:zE);
    %closing (dilate then erode) fills gaps and smooths the surface
    vesselCrop = imclose(vesselCrop,strel('sphere',I-1));
    %opening (erode then dilate) removes corner protrusions and small artifacts
    vesselCrop = imopen( vesselCrop,strel('sphere',1  ));
    %write back to non-cropped nifti
    mrixVesselSeg{v}.vol(xS:xE, yS:yE, zS:zE) = vesselCrop;
    MRIwrite(mrixVesselSeg{v}, mrixVesselSeg{v}.fspec);
    
    %get surf
    vesselSurf{v} = replace(mrixVesselSeg{v}.fspec, '.nii', '.vtk');
    vmtk_surfFromSeg(mrixVesselSeg{v}.fspec, vesselSurf{v})

    vesselCenterlines{v} = replace(vesselSurf{v}, '/seg/', '/surf/');
    cmd = {src.vmtk};
    cmd{end+1} = 'vmtkcenterlines \';
    cmd{end+1} = ['-ifile ' vesselSurf{v} ' \'];
    cmd{end+1} = ['-ofile ' vesselCenterlines{v}];
    disp(strjoin(cmd,newline));

    % vmtk_surfSmoothing(vesselSurf{v}, vesselSurf{v}, options)
    vmtk_viewVolAndSurf(mrixSeg.fspec, vesselSurf{v})
    vmtk_viewVolAndSurf(mrixSeg.fspec, vesselCenterlines{v})
    










    mriVesselSeg{v} = MRIread(mriVesselSeg{v}.fspec);

    
    mriVesselSegCropped{v}.fspec = replace(mriVesselSeg{v}.fspec, '.nii.gz', ['_cropped.nii.gz']);
    cmd = {};
    cmd{end+1} = ['python ' fullfile(toolDir, 'bassWrap-reg', 'crop_nifti_ants.py') ' \'];
    cmd{end+1} = [mriVesselSeg{v}.fspec ' \'];
    cmd{end+1} = [mriVesselSegCropped{v}.fspec];
    system(strjoin(cmd,newline),'-echo');

    mriVesselSeg{v}        = MRIread(mriVesselSeg{v}.fspec);
    mriVesselSegCropped{v} = MRIread(mriVesselSegCropped{v}.fspec);
    size(mriVesselSeg{v}.vol)
    size(mriVesselSegCropped{v}.vol)

    freeview /scratch/users/Proulx-S/vesselReg3/tmp/tof/seg/tof_vessel006.nii.gz    /scratch/users/Proulx-S/vesselReg3/tmp/tof/seg/tof_vessel006_cropped.nii.gz

    voxSz1 = [mriVesselSeg{v}.xsize mriVesselSeg{v}.ysize mriVesselSeg{v}.zsize]';
    voxSz2 = [mriVesselSegCropped{v}.xsize mriVesselSegCropped{v}.ysize mriVesselSegCropped{v}.zsize]';
    mriVesselSeg{v}.tkrvox2ras(1:end-1,end)./voxSz1
    ceil(mriVesselSegCropped{v}.tkrvox2ras(1:end-1,end)./voxSz2)



    [~, bbox] = auto_crop_nifti(mriVesselSeg{v}.fspec, mriVesselSegCropped{v}.fspec,[],1);
    mriVesselSegCropped{v} = MRIread(mriVesselSegCropped{v}.fspec);

    mriVesselSegCropped{v}.tkrvox2ras
    x = bbox.min(1):bbox.max(1);
    y = bbox.min(2):bbox.max(2);
    z = bbox.min(3):bbox.max(3);
    tmp = mriVesselSeg{v}.vol(y+1, x+1, z+1)-mriVesselSegCropped{v}.vol;

    

    mriVesselSegCroppedUpsampled{v}.fspec = replace(mriVesselSeg{v}.fspec, '.nii.gz', ['_cropped_upsampled.nii.gz']);
    


    mriVesselSegUpSampled{v}.fspec = replace(mriVesselSeg{v}.fspec, '.nii.gz', ['_upsampled' num2str(I) '.nii.gz']);
    cmd = {src.fs};
    cmd{end+1} = 'mri_convert \';
    cmd{end+1} = ['--resample_type nearest \'];
    cmd{end+1} = ['--upsample ' num2str(I) ' \'];
    cmd{end+1} = [mriVesselSeg{v}.fspec ' \'];
    cmd{end+1} = [mriVesselSegUpSampled{v}.fspec];
    % system(strjoin(cmd,newline),'-echo');
    % cmd = {src.fs};
    cmd{end+1} = 'mri_morphology \';
    cmd{end+1} = [mriVesselSegUpSampled{v}.fspec ' \'];
    cmd{end+1} = ['dilate ' num2str(I-1) ' \'];
    cmd{end+1} = [replace(mriVesselSegUpSampled{v}.fspec, '.nii.gz', '_dilate.nii.gz')];
    cmd{end+1} = 'mri_morphology \';
    cmd{end+1} = [mriVesselSegUpSampled{v}.fspec ' \'];
    cmd{end+1} = ['erode ' num2str(I-1) ' \'];
    cmd{end+1} = [replace(mriVesselSegUpSampled{v}.fspec, '.nii.gz', '_erode.nii.gz')];
    % system(strjoin(cmd,newline),'-echo');
    disp(strjoin(cmd,newline));

    mriVesselSeg{v} = MRIread(mriVesselSeg{v}.fspec);
    mriVesselSegUpSampled{v} = MRIread(mriVesselSegUpSampled{v}.fspec);
    mriVesselSegUpSampledClosed{v} = MRIread(replace(mriVesselSegUpSampled{v}.fspec, '.nii.gz', '_close.nii.gz'));

    figure; ax = {};
    ht = tiledlayout(3,3); ht.TileSpacing = 'compact'; ht.Padding = 'compact';
    for i = 1:3
    ax{end+1} = nexttile;
    imagesc([1 400],[1 400],mriVesselSeg{v}.vol(:,:,end/2+i-1)); axis image off square
    ax{end+1} = nexttile;
    imagesc([1 400],[1 400],mriVesselSegUpSampled{v}.vol(:,:,end/2+i-1)); axis image off square
    ax{end+1} = nexttile;
    imagesc([1 400],[1 400],mriVesselSegUpSampledClosed{v}.vol(:,:,end/2+i-1)); axis image off square
    end
    drawnow;
    set([ax{:}],'YLim',[100 175],'XLim',[50 150]+20);
    colormap gray

    cmd = {src.fs};
    cmd{end+1} = 'freeview \';
    cmd{end+1} = ['-v ' mriVesselSegUpSampled{v}.fspec ' \'];
    cmd{end+1} = ['-v ' replace(mriVesselSegUpSampled{v}.fspec, '.nii.gz', '_close.nii.gz') ' \'];
    cmd{end+1} = ['-v ' mriVesselSeg{v}.fspec ':resample=nearest'];
    disp(strjoin(cmd,newline));



    mriVesselSegUpSampled{v} = MRIread(mriVesselSegUpSampled{v}.fspec);
    for i = 1:I
        mriVesselSegUpSampled{v}.vol = imerode(  imdilate(mriVesselSegUpSampled{v}.vol,strel('sphere',1))  ,strel('sphere',1));
    end
    mriVesselSegUpSampled{v}.fspec = replace(mriVesselSeg{v}.fspec, '.nii.gz', ['_upsampled' num2str(I) '.nii.gz']);
    MRIwrite(mriVesselSegUpSampled{v}, mriVesselSegUpSampled{v}.fspec);    
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Skeletonize using scikit-image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use scikit-image's skeletonize (Lee method for 3D)
tof_skeleton_skimage = fullfile(fileparts(tof_vesselSeg), 'tof_skeleton_skimage.nii.gz');
skeletonize_nifti(tof_vesselSeg, tof_skeleton_skimage, 'lee', forceThis);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mri = MRIread(tof_skeleton_skimage);
% CC = bwconncomp(mri.vol,26);








tof_skeleton_skimage_label = fullfile(fileparts(tof_skeleton_skimage), 'tof_skeleton_label.nii.gz');
if forceThis || ~exist(tof_skeleton_skimage_label,'file')
    label_connected_components_nifti(tof_skeleton_skimage, tof_skeleton_skimage_label, 3);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



forceThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write mask of selected vessels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CC = bwconncomp(mriSegAll.vol,26);

cmd = {src.fs};
cmd{end+1} = 'freeview \';
cmd{end+1} = ['-v ' tof ' \'];
cmd{end+1} = ['-v ' tof_skeleton_skimage ' \'];
cmd{end+1} = ['-v ' tof_skeleton_skimage_label];
disp(strjoin(cmd,newline));

okCompIdx = [31 20 18 25 15 19 8 56 14 29 50 6 9];
mriVesselLabels = MRIread(tof_skeleton_skimage_label);
vesselMaskList = cell(length(okCompIdx), 1);
for i = 1:length(okCompIdx)
    vesselMaskList{i} = fullfile(fileparts(tof_skeleton_skimage_label), ['tof_skeleton_label_' num2str(okCompIdx(i)) '.nii.gz']);
    if exist(vesselMaskList{i}, 'file') && ~forceThis; continue; end
    tmp = mriVesselLabels;
    tmp.vol = zeros(size(tmp.vol));
    tmp.vol(mriVesselLabels.vol == okCompIdx(i)) = 1;
    MRIwrite(tmp, vesselMaskList{i});



    mriSkel = tmp;
    mriSegAll = MRIread(tof_vesselSeg);
    mriVessel = mriSkel;
    % Extract vessel from full segmentation using skeleton
    mriVessel.vol = double(extract_vessel_from_segmentation(mriSegAll.vol, mriSkel.vol));
    mriVessel.vol(mriSkel.vol==1) = 0.5;
    MRIwrite(mriVessel, replace(vesselMaskList{i},'.nii.gz','_fullVessel.nii.gz'));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Upsample vessel mask and dilate





return

forceThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract centerlines from skeleton mask of each vessel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create centerline VTK files for each identified vessel skeleton
% Using graph extraction to preserve full branching structure
% Use original TOF volume as reference for affine transformation to ensure
% proper spatial alignment
vesselCenterlineList = cell(size(vesselMaskList));
for i = 1:length(vesselMaskList)
    % Extract graph from skeleton component and write VTK centerline
    % This creates a full graph where each skeleton voxel is a node
    % and edges connect adjacent voxels, preserving branching structure
    vesselCenterlineList{i} = replace(replace(vesselMaskList{i}, '.nii.gz', '.vtk'),'/seg/','/centerlines/');
    if exist(vesselCenterlineList{i}, 'file') && ~forceThis; continue; end
    skeleton_to_graph_vtk(vesselMaskList{i}, vesselCenterlineList{i}, 26, forceThis);
    % skan_skeletonMask_to_vtk(vesselMaskList{i}, vesselCenterlineList{i}, forceThis);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% i = 3;
vmtk_viewVolAndSurf(vesselMaskList{i}, vesselCenterlineList{i});






%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit vessel tube model
%%%%%%%%%%%%%%%%%%%%%%%%
vesselRefinedCenterlineList{i} = replace(vesselCenterlineList{i}, '.vtk', '_refined.vtk');

cmd        = {src.vmtk};
cmd{end+1} = 'vmtkactivetubes \';
cmd{end+1} = ['-imagefile ' tof ' \'];
cmd{end+1} = ['-ifile ' vesselCenterlineList{i} ' \'];
cmd{end+1} = ['-ofile ' vesselRefinedCenterlineList{i}];
system(strjoin(cmd,newline),'-echo');

vmtk_viewVolAndSurf(tof, vesselCenterlineList{i});
vmtk_viewVolAndSurf(tof, vesselRefinedCenterlineList{i});



forceThis = 1;
vesselRefinedCenterlineList = cell(size(vesselCenterlineList));
for i = 1:length(vesselCenterlineList)
    % Refine centerlines using vmtkactivetubes to fit vessel tube model
    vesselRefinedCenterlineList{i} = replace(vesselCenterlineList{i}, '.vtk', '_refined.vtk');
    if exist(vesselRefinedCenterlineList{i}, 'file') && ~forceThis; continue; end
    options = struct();
    options.iterations = 100;
    options.potentialweight = 1.0;
    options.stiffnessweight = 1.0;
    options.forceThis = forceThis;
    vmtk_activetubes(tof, vesselCenterlineList{i}, vesselRefinedCenterlineList{i}, options);
end


%% %%%%%%%%%%%%%%%%%%%%%

