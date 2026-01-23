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

    S=2;
    %tof
    in = fullfile(roi{S}.(acq).(task).rCond.volAnat.tof.folder,roi{S}.(acq).(task).rCond.volAnat.tof.name);
    out = fullfile(projectCode,'data','tof'); if ~exist(out,'dir'); mkdir(out); end
    tof = fullfile(out,'tof.nii.gz');
    copyfile(in,tof);

    %vfMRI
    in = roi{S}.(acq).(task).rCond.fOrigList{1};
    out = fullfile(projectCode,'data','vfMRI'); if ~exist(out,'dir'); mkdir(out); end
    vfMRI = fullfile(out,'vfMRI.nii.gz');
    copyfile(in,vfMRI);
    %% %%%%%%%%%%%%%%%
else
    %%%%%%%%%%%%%%%%%
    %% Get data files
    tof   = fullfile(projectCode,'data','tof'  ,'tof.nii.gz'  );
    vfMRI = fullfile(projectCode,'data','vfMRI','vfMRI.nii.gz');
    %% %%%%%%%%%%%%%%
end


forceThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vesselboost vessel segmentation of tof volume data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tof_vesselSeg = fullfile(projectScratch,'tof','seg'); if ~exist(tof_vesselSeg,'dir'); mkdir(tof_vesselSeg); end
tof_vesselSeg = fullfile(tof_vesselSeg,'tof.nii.gz');
if forceThis || ~exist(tof_vesselSeg,'file')
    vesselboost_prediction(fileparts(tof),fileparts(tof_vesselSeg),vesselBoostModel,4);
end
tof_vesselSegX = replace(tof_vesselSeg, '.nii.gz', 'US.nii');
%reduce precision and don't compress
cmd = {src.fs};
cmd{end+1} = 'mri_convert \';
cmd{end+1} = '--out_data_type uchar \';
cmd{end+1} = [tof_vesselSeg ' \'];
cmd{end+1} = tof_vesselSegX;
system(strjoin(cmd,newline),'-echo');
%upsample
I = 8;
cmd = {src.fs};
cmd{end+1} = 'mri_convert \';
cmd{end+1} = ['--upsample ' num2str(I) ' \'];
cmd{end+1} = ['--resample_type nearest \'];
cmd{end+1} = [tof_vesselSegX ' \'];
cmd{end+1} = tof_vesselSegX;
system(strjoin(cmd,newline),'-echo');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Individualize and skeletonize vessels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find connected components (vessels)
mrixSeg = MRIread(tof_vesselSegX);
mrixSeg.vol = logical(mrixSeg.vol);
disp('bwconncomp: running...');
CC = bwconncomp(mrixSeg.vol,26);
disp('bwconncomp: done');
% create a mask for each vessel
mrixVesselSeg = cell(CC.NumObjects, 1);
mrixVesselSkel = cell(CC.NumObjects, 1);
n = zeros(CC.NumObjects, 1);
for v = 1:CC.NumObjects
    %!!!!should crop here!!!!!
    mrixVesselSeg{v} = mrixSeg;
    mrixVesselSeg{v}.vol = false(size(mrixSeg.vol));
    mrixVesselSeg{v}.vol(CC.PixelIdxList{v}) = true;
    % skeletonize individual vessel
    mrixVesselSkel{v} = mrixSeg;
    disp(['bwskel: running... ' num2str(v) '/' num2str(CC.NumObjects)]);
    mrixVesselSkel{v}.vol = bwskel(mrixVesselSeg{v}.vol);
    nx(v) = nnz(mrixVesselSkel{v}.vol);
end
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

