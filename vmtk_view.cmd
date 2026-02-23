ml vmtk/1.5.0
vmtkimagereader -ifile /scratch/users/Proulx-S/vesselReg3/msVesselBoost/tof/scale-1_tof.nii.gz \
--pipe vmtksurfacereader -ifile /scratch/users/Proulx-S/vesselReg3/vesselPrc/centerlines/vessel-00_centerline.vtk \
--pipe vmtkrenderer \
--pipe vmtkimageviewer -i @vmtkimagereader.o -display 0 -textureinterpolation 0 \
--pipe vmtksurfaceviewer -color 0 1 0 -display 0 \
--pipe vmtksurfacereader -ifile /scratch/users/Proulx-S/vesselReg3/vesselPrc/centerlines/vessel-00_desc-activetubes_centerline.vtk \
--pipe vmtksurfaceviewer -i @vmtksurfacereader.o -color 1 0 0 -display 1
