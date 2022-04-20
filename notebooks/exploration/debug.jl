module TestPatterns
using HDF5, ReusePatterns
# import ReusePatterns: @forward

struct DsetWrapper
    dset::HDF5.Dataset
end

@forward((DsetWrapper,:dset),HDF5.Dataset)
end