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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get msVesselBoost pointers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msVesselBoost_pointerFile = fullfile(info.project.code,'msVesselBoost',['sub-' num2str(S) '_vesselBoostPointer.mat']);
load(msVesselBoost_pointerFile); % overwrites info, but should be fine
info.pointerFile.msVesselBoost = info.pointerFile.vessel;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%
info;




forceThis = 1;
consensusThreshValue = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Label individualizable vessels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fMask  = fullfile(info.project.code, 'tmp', 'mask' , ['scale-C_thresh-' num2str(consensusThreshValue) '_mask.nii.gz']);
fLabel = fullfile(info.project.code, 'selectVessels', 'label',  'scale-C_label.nii.gz');
if ~exist(fileparts(fMask),'dir'); mkdir(fileparts(fMask)); end
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
else
    disp('thresholding and splitting vessels... already done');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fMask;
fLabel;
info.subject.label.fList = {fLabel};





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare vfMRI data if available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vfMRI_in = dir(fullfile(info.project.code,'data','vfMRI','vfMRI.nii.gz'));
if ~isempty(vfMRI_in)
    % Average the vfMRI timeseries across time using ants' AverageImages
    vfMRI_in  = fullfile(vfMRI_in.folder, vfMRI_in.name);
    vfMRI_out = fullfile(info.project.code, 'selectVessels', 'vfMRI', 'vfMRI.nii.gz');
    if ~exist(fileparts(vfMRI_out),'dir'); mkdir(fileparts(vfMRI_out)); end
    if ~exist(vfMRI_out, 'file')
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



return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% View vessel label map in freeview and select vessel labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmd = {src.fs};
cmd{end+1} = 'freeview \';
valRange = [0 length(info.subject.seg.fList)-1]; valRange(valRange==0) = valRange(valRange==0) + 0.01;
cmd{end+1} = ['-v ' info.subject.label.fList{1} ':colormap=heat:heatscale=' num2str(valRange(1)) ',' num2str(valRange(2)) ' \']; % concensus segmentation
cmd{end+1} = ['-v ' info.subject.tof.fList{1} ':resample=nearest \']; % original resolution tof
cmd{end+1} = ['-v ' info.subject.vfMRI.fList{1}]; % vfMRI

info.subject.visualize(end+1).f = fullfile(fileparts(fileparts(info.subject.label.fList{1})),'selectVessels_view.cmd');
if forceThis || ~exist(info.subject.visualize(end).f,'file')
    fid = fopen(info.subject.visualize(end).f,'w');
    fprintf(fid,'%s\n',strjoin(cmd,newline));
    fclose(fid);
end
info.subject.visualize(end).cmd = cmd';
%% %%%%%%%%%%%%%%%%%
info.subject.visualize(end)
disp(strjoin(replace(info.subject.visualize(end).cmd,'/scratch/users/Proulx-S','/scratch/Proulx-S'),newline));











forceThis = 0;
% Note: we might want to consider a more strict criteria (connectivity parameter for bwconncomp) to better separate vessels that are close to one another.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Individualize vessels into a single-vessel label map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vesselIdx = sort([7 9 47 19 27 23 17 15 8 16 100 31 28 12 13 48 26 11 20 36 14]);

fMask  = fullfile(info.project.code, 'tmp', 'mask' , ['scale-C_thresh-' num2str(threshValue) '_mask.nii.gz']);
fLabel = fullfile(info.project.code, 'tmp', 'label',  'scale-C_label.nii.gz');
if ~exist(fileparts(fMask),'dir'); mkdir(fileparts(fMask)); end
if ~exist(fileparts(fLabel),'dir'); mkdir(fileparts(fLabel)); end
if forceThis || ~exist(fMask,'file') || ~exist(fLabel,'file')
    % Create vessel mask by thresholding consensus segmentation)
    disp('thresholding consensus segmentation... computing');
    if ~exist('mriSum','var'); mriSum = MRIread(tofCropPrcSclSegUpsSum); end
    mriMask = mriSum;
    mriMask.fspec = fMask;
    mriMask.vol = mriSum.vol >= threshValue;
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

    % % visualize vessel labels and select vessels (vesselIdx)
    % MRIwrite(mriLabel, fLabel,'ushort');
    % cmd = {src.fs};
    % cmd{end+1} = 'freeview \';
    % cmd{end+1} = ['-v ' fLabel                            ' \'];
    % cmd{end+1} = ['-v ' vfMRI            ':resample=nearest \'];
    % cmd{end+1} = ['-v ' tofCropPrcScl{1} ':resample=nearest'  ];
    % disp(strjoin(cmd,newline));
    
    % write mask and label files with only selected vessels
    mriLabel.vol(~ismember(mriLabel.vol,vesselIdx)) = 0;
    mriMask.vol( ~ismember(mriMask.vol ,vesselIdx)) = 0;
    disp('saving...');
    MRIwrite(mriLabel, fLabel,'ushort');
    MRIwrite(mriMask , fMask , 'uchar');
    disp('thresholding and splitting vessels... done');
else
    disp('thresholding and splitting vessels... already done');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fMask;
fLabel;











forceThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Crop out each single vessel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
nVessel  = length(vesselIdx);
fTofList = tofCropPrcScl([1 end]);
fSegList = [tofCropPrcSclSegUps {tofCropPrcSclSegUpsSum}];
fMask;
fLabel;
% outputs
fVesselTofList   = cell(nVessel, length(fTofList));
fVesselSegList   = cell(nVessel, length(fSegList));
fVesselMaskList  = cell(nVessel, 1);
fVesselLabelList = cell(nVessel, 1);
for v = 1:nVessel
    fVesselTofList{v,1} = fullfile(info.project.code, 'result', 'tof', ['vessel-' sprintf('%02d',v) '_scale-' num2str(min(scaleList)) '_tof.nii.gz']);
    fVesselTofList{v,2} = fullfile(info.project.code, 'result', 'tof', ['vessel-' sprintf('%02d',v) '_scale-' num2str(max(scaleList)) '_tof.nii.gz']);
    for s = 1:length(scaleList)
        fVesselSegList{v,s} = fullfile(info.project.code, 'result', 'seg', ['vessel-' sprintf('%02d',v) '_scale-' num2str(scaleList(s)) '_seg.nii.gz']);
    end
    fVesselSegList{v,end} = fullfile(info.project.code, 'result', 'seg'  , ['vessel-' sprintf('%02d',v) '_scale-C_seg.nii.gz']);
    fVesselMaskList{v}    = fullfile(info.project.code, 'result', 'mask' , ['vessel-' sprintf('%02d',v) '_scale-C_mask.nii.gz']);
    fVesselLabelList{v}   = fullfile(info.project.code, 'result', 'label', ['vessel-' sprintf('%02d',v) '_scale-C_label.nii.gz']);
end

if forceThis || any(~cellfun(@(f) exist(f, 'file'), [fVesselTofList(:); fVesselSegList(:); fVesselMaskList(:); fVesselLabelList(:)]))
    
    if ~exist('mriLabel','var'); mriLabel = MRIread(fLabel); end
    mriMask_scaleC = mriLabel;
    mriMask_scale1 = MRIread(tofCropPrcScl{1},1);
    
    for vIdx = 1:nVessel
        v = vesselIdx(vIdx);

        % Write temporary vessel masks
        mriMask_scale1.fspec = fullfile(info.project.code, 'tmp', 'mask', ['vessel-' sprintf('%02d',v) '_scale-1_tmpMask.nii.gz']);
        mriMask_scaleC.fspec = fullfile(info.project.code, 'tmp', 'mask', ['vessel-' sprintf('%02d',v) '_scale-C_tmpMask.nii.gz']);
        mriMask_scaleC.vol = mriLabel.vol==v;
        mriMask_scale1.vol = imresize3(mriMask_scaleC.vol, mriMask_scale1.volsize, 'box'    );
        mriMask_scaleC.vol = imresize3(mriMask_scale1.vol, mriMask_scaleC.volsize, 'nearest');

        MRIwrite(mriMask_scale1, mriMask_scale1.fspec, 'uchar');
        MRIwrite(mriMask_scaleC, mriMask_scaleC.fspec, 'uchar');

        % Crop vessels
        % tof
        for s = 1:length(fTofList)
            if ~exist(fileparts(fVesselTofList{v,s}),'dir'); mkdir(fileparts(fVesselTofList{v,s})); end
            if forceThis || ~exist(fVesselTofList{v,s},'file')
                disp(['cropping vessel (' num2str(vIdx) '/' num2str(nVessel) ') from tof (' num2str(s) '/' num2str(length(fTofList)) ')... computing']);
                cmd = {src.ants};
                cmd{end+1} = 'ExtractRegionFromImageByMask 3 \';
                cmd{end+1} = [fTofList{s} ' \'];
                cmd{end+1} = [fVesselTofList{v,s} ' \'];
                scaleStr = strsplit(fTofList{s},filesep); scaleStr = scaleStr{end}; scaleStr = strsplit(scaleStr,'_'); scaleStr = scaleStr{contains(scaleStr,'scale-')};
                switch scaleStr
                    case 'scale-1'
                        cmd{end+1} = [mriMask_scale1.fspec ' \'];
                        cmd{end+1} = ['1 1']; % [label, padRadius]
                    case {'scale-8','scale-C'}
                        cmd{end+1} = [mriMask_scaleC.fspec ' \'];
                        cmd{end+1} = ['1 8']; % [label, padRadius]
                    otherwise
                        error('Invalid scale');
                end
                [status,result] = system(strjoin(cmd, newline), '-echo');
                disp(['cropping vessel (' num2str(vIdx) '/' num2str(nVessel) ') from tof (' num2str(s) '/' num2str(length(fTofList)) ')... done']);
            else
                disp(['cropping vessel (' num2str(vIdx) '/' num2str(nVessel) ') from tof (' num2str(s) '/' num2str(length(fTofList)) ')... already done']);
            end
        end

        % seg
        for s = 1:length(fSegList)
            if ~exist(fileparts(fVesselSegList{v,s}),'dir'); mkdir(fileparts(fVesselSegList{v,s})); end
            if forceThis || ~exist(fVesselSegList{v,s},'file')
                disp(['cropping vessel (' num2str(vIdx) '/' num2str(nVessel) ') from seg (' num2str(s) '/' num2str(length(fSegList)) ')... computing']);
                cmd = {src.ants};
                cmd{end+1} = 'ExtractRegionFromImageByMask 3 \';
                cmd{end+1} = [fSegList{s} ' \'];
                cmd{end+1} = [fVesselSegList{v,s} ' \'];
                cmd{end+1} = [mriMask_scaleC.fspec ' \'];
                cmd{end+1} = ['1 8']; % [label, padRadius]
                [status,result] = system(strjoin(cmd, newline), '-echo');
                disp(['cropping vessel (' num2str(vIdx) '/' num2str(nVessel) ') from seg (' num2str(s) '/' num2str(length(fSegList)) ')... done']);
            else
                disp(['cropping vessel (' num2str(vIdx) '/' num2str(nVessel) ') from seg (' num2str(s) '/' num2str(length(fSegList)) ')... already done']);
            end
        end
        % mask
        if ~exist(fileparts(fVesselMaskList{v}),'dir'); mkdir(fileparts(fVesselMaskList{v})); end
        if forceThis || ~exist(fVesselMaskList{v},'file')
            disp(['cropping vessel (' num2str(vIdx) '/' num2str(nVessel) ') from mask... computing']);
            cmd = {src.ants};
            cmd{end+1} = 'ExtractRegionFromImageByMask 3 \';
            cmd{end+1} = [fMask ' \'];
            cmd{end+1} = [fVesselMaskList{v} ' \'];
            cmd{end+1} = [mriMask_scaleC.fspec ' \'];
            cmd{end+1} = ['1 8']; % [label, padRadius]
            [status,result] = system(strjoin(cmd, newline), '-echo');
            disp(['cropping vessel (' num2str(vIdx) '/' num2str(nVessel) ') from mask... done']);
        else
            disp(['cropping vessel (' num2str(vIdx) '/' num2str(nVessel) ') from mask... already done']);
        end
        %label
        if ~exist(fileparts(fVesselLabelList{v}),'dir'); mkdir(fileparts(fVesselLabelList{v})); end
        if forceThis || ~exist(fVesselLabelList{v},'file')
            disp(['cropping vessel (' num2str(vIdx) '/' num2str(nVessel) ') from label... computing']);
            cmd = {src.ants};
            cmd{end+1} = 'ExtractRegionFromImageByMask 3 \';
            cmd{end+1} = [fLabel ' \'];
            cmd{end+1} = [fVesselLabelList{v} ' \'];
            cmd{end+1} = [mriMask_scaleC.fspec ' \'];
            cmd{end+1} = ['1 8']; % [label, padRadius]
            [status,result] = system(strjoin(cmd, newline), '-echo');
            disp(['cropping vessel (' num2str(vIdx) '/' num2str(nVessel) ') from label... done']);
        else
            disp(['cropping vessel (' num2str(vIdx) '/' num2str(nVessel) ') from label... already done']);
        end
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
fVesselTofList;
fVesselSegList;
fVesselMaskList;
fVesselLabelList;

return

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