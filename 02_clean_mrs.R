library(spant)
library(RNifti)
library(fslr)
library(fs)

# get the directory structure of the raw data
dirs <- dir_ls(recurse = 2, glob = "MRS_MRI_CLEANED/*/*", type = "directory")

# copy the input directory structure and make a new one for storing results
out_dirs <- path("MRS_MRI_PROC", path_rel(dirs, "MRS_MRI_CLEANED"))
dir_create(out_dirs)

# find the mrs raw data files and copy to the new directory structure
mrs_files <- Sys.glob("MRS_MRI_CLEANED/*/*/*.dat")
mrs_file_size <- file.info(mrs_files)$size

# this distinguishes between ws and wref
ws <- mrs_file_size > 2000000
wref <- !ws

# copy the ws files
ws_dest <- path(out_dirs, "ws.dat")
file_copy(mrs_files[ws], ws_dest)
wref_dest <- path(out_dirs, "wref.dat")
file_copy(mrs_files[wref], wref_dest)