library(spant)
library(RNifti)
library(fs)

# find the data directories
data_dirs <- dir_ls(glob = "MRS_MRI_PROC/*/*", recurse = 2)

mri_path             <- path(data_dirs, "t1.nii.gz")
seg_path             <- path(data_dirs, "t1_seg.nii.gz")
ws_path              <- path(data_dirs, "ws.dat")
wref_path            <- path(data_dirs, "wref.dat")
vox_overlay_path     <- path(data_dirs, "vox_overlay.png")
vox_overlay_seg_path <- path(data_dirs, "vox_overlay_seg.png")

# generate vox pos files
mapply(plot_voi_overlay, mri = mri_path, voi = ws_path,
       export_path = vox_overlay_path)

# generate vox PVC files and data
pvs <- mapply(plot_voi_overlay_seg, mri_seg = seg_path, voi = ws_path, 
              export_path = vox_overlay_seg_path)

# write the partial volume numbers
write.table(pvs, "pvs.txt")

# write the data directories to keep later processing steps consistent
write.table(data_dirs, "data_dirs.txt")