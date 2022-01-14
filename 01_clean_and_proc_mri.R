library(spant)
library(RNifti)
library(fslr)
library(fs)

# get the directory structure of the raw data
dirs <- dir_ls(recurse = 2, glob = "MRS_MRI_CLEANED/*/*", type = "directory")

# copy the input directory structure and make a new one for storing results
out_dirs <- path("MRS_MRI_PROC", path_rel(dirs, "MRS_MRI_CLEANED"))
dir_create(out_dirs)

# find the T1 files and copy to the new directory structure
t1_files <- dir_ls(glob = "MRS_MRI_CLEANED/*/*/*nii.gz", recurse = 3)
t1_dest <- path(out_dirs, "t1.nii.gz")
file_copy(t1_files, t1_dest)

# run bet on all files
brain_dest <- path(out_dirs, "brain")
mapply(fslbet, infile = t1_dest, outfile = brain_dest,
       MoreArgs = list(retimg = FALSE, opts = "-R", betcmd = "bet"))

# run fast on all files
fast_dest <- path(out_dirs, "t1")
mapply(fslfast, file = brain_dest, outfile = fast_dest,
       MoreArgs = list(retimg = FALSE))