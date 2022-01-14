library(spant)
library(RNifti)
library(fs)

# find the data directories
data_dirs <- dir_ls(glob = "MRS_MRI_PROC/*/*", recurse = 2)

mri_path             <- path(data_dirs, "t1.nii.gz")
seg_path             <- path(data_dirs, "t1_seg.nii.gz")
ws_path              <- path(data_dirs, "ws.dat")
w_ref_path           <- path(data_dirs, "wref.dat")
vox_overlay_path     <- path(data_dirs, "vox_overlay.png")
vox_overlay_seg_path <- path(data_dirs, "vox_overlay_seg.png")

temp  <- downsample_mrs_td(read_mrs(ws_path[1]))

# basis simulation
basis <- sim_basis_1h_brain_press(temp, TE1 = 0.01, TE2 = 0.07)

idx <- seq_len(length(ws_path))

opts <- abfit_opts()

extra <- data.frame(scan_id = substr(ws_path, 25, 34),
                    subj    = substr(ws_path, 22, 23))

res_full <- svs_1h_brain_batch_analysis(ws_path[idx], w_ref_path[idx],
                                  mri_seg = seg_path[idx], extra = extra[idx,],
                                  basis = basis, te = 0.08, tr = 3)

saveRDS(res_full, "all_fits_simple_basis.rds")