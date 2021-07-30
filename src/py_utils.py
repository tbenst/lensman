import cv2
import numpy as np
import tables
from tqdm import tqdm

def assign_tseries_to_dset(dset, tseries, channel=0, batch_size=1024):
    "dset is 5D (T, C, Z, H, W), tseries is 4D (T, Z, H, W)"
    if len(tseries.shape) == 4:
        T, Z, H, W = tseries.shape
        for z in range(Z):
            print('Saving z-plane {} out of {}'.format(z + 1, Z))
            # Using batching to write out tseries (these blobs are massive ~75GB)
            for b in tqdm(range(0, T, batch_size)):
                dset[b:b + batch_size, channel, z] = tseries[b:b + batch_size, z]
    else:
        raise NotImplementedError


def resize_images(images, target_width, target_height, interpolation=cv2.INTER_LINEAR):
    assert len(images.shape) == 4
    _, _, original_height, original_width = images.shape
    fx = target_width / original_width
    fy = target_height / original_height

    is_binary_image = (images.dtype == np.uint8)
    if is_binary_image:
        assert interpolation == cv2.INTER_NEAREST

    example_resized_image = cv2.resize((images[0, 0]), None, fx=fx, fy=fy, interpolation=interpolation)
    assert example_resized_image.shape == (target_height, target_width), 'Only allowing exact resizes'

    resized_images = np.zeros((images.shape[0], images.shape[1], target_height, target_width), dtype=np.float32)
    for t, volume in enumerate(images):
        for z, image in enumerate(volume):
            resized_images[t, z] = cv2.resize(image, None, fx=fx, fy=fy, interpolation=interpolation)

    # If images passed in were boolean then return as such
    if is_binary_image:
        assert set(np.unique(resized_images)).issubset({0, 1})
        resized_images = np.uint8(resized_images)

    return resized_images


def write_experiment_to_tyh5(
    tseries,
    output_path,
    stim_masks=None,
    stim_used_at_each_timestep=None,
    num_channels=1,
    compression_level=3,
    small_size=256
):
    """Write out a full experiment (including stim_masks if available) to .ty.h5

    Structure written out is:

    IMAGING
    |--- RAW
    |--- SMALL
    STIM
    |--- MASKS
        |--- RAW
        |--- SMALL
    |--- STIM_USED_AT_EACH_TIMESTEP


    Parameters
    ----------
    tseries: np.uint16 array
        (H, W, Z, T) time series

    output_path: str
        Path to write .ty.h5 to

    stim_masks: np.bool array, default: None
        Set of unique binary masks used for stimulation
        Shape: (n_stimuli, Z, H, W)

    stim_used_at_each_timestep: np.int64 array, default: None
        Entry i is 0 if no stimulation was active at ith timestep, otherwise
        stim used at this timestep is stim_masks[i - 1]
        Shape: (T,)

    num_channels: int, default: 1
        This isn't really doing much as we currently do not have support for more channels
        TODO(allan.raventos): add support for multiple channels

    compression_level: int, default: 3
        Compression level to use for 'blosc:zstd' when writing out data

    small_size: int, default: 256
        Size to write to '/imaging/small', which will generally be used for training

    Notes
    -----
    This function is less useful if full tseries doesn't fit in RAM (in that case we'll
    need to load from TIFF files directly - can write this later if/when needed)

    `nothing` on Julia side converts to `None` here as desired
    """
    if tseries.dtype != np.uint16 or len(tseries.shape) != 4:
        raise TypeError('tseries does not have the expected type or shape')

    H, W, Z, T = tseries.shape

    # Check that `stim_masks` has the correct type and shape
    if stim_masks is not None:
        if stim_masks.dtype != np.bool:
            raise TypeError('stim_masks does not have the expected type')
        if len(stim_masks.shape) != 4 or stim_masks.shape[1:] != (Z, H, W):
            raise ValueError('stim_masks does not have the expected dimensions')

        # Hacky but seemingly necessary, h5py won't read the default pytables bool type
        stim_masks = np.uint8(stim_masks)
        assert set(np.unique(stim_masks)).issubset({0, 1})

        # Check that `stim_used_at_each_timestep` has the correct type and shape
        if stim_used_at_each_timestep is not None and stim_used_at_each_timestep.dtype != np.int64:
            raise TypeError('stim_used_at_each_timestep does not have the expected type')
        if stim_used_at_each_timestep.shape != (T, ):
            raise ValueError('stim_used_at_each_timestep should be as long as tseries')

    with tables.open_file(output_path, 'w') as tyh5:

        compression_filter = tables.filters.Filters(complevel=compression_level, complib='blosc:zstd')

        tseries = np.transpose(tseries, (3, 2, 0, 1))  # HWZT -> TZHW
        imaging_group = tyh5.create_group('/', 'imaging')
        imaging_attrs = tyh5.root['imaging']._v_attrs
        imaging_attrs['num_frames'] = T
        imaging_attrs['num_z_planes'] = Z

        # Create raw resolution
        dset_shape = (T, num_channels, Z, H, W)
        dset_raw = tyh5.create_carray(
            imaging_group, 'raw', tables.Atom.from_dtype(tseries.dtype), shape=dset_shape, filters=compression_filter
        )
        print('Writing tseries (TZHW={}) to "/imaging/raw" (TCZHW={})'.format(tseries.shape, dset_raw.shape))
        assign_tseries_to_dset(dset_raw, tseries)
        dset_raw.attrs['dimensions'] = 'TCZHW'

        # Generate '/imaging/small' at the requested size
        # TB comment: perhaps convenient if small is always 256, regardless of raw
        dset_small_shape = (T, num_channels, Z, small_size, small_size)
        tseries_small = resize_images(tseries, small_size, small_size, interpolation=cv2.INTER_LINEAR)
        dset_small = tyh5.create_carray(
            imaging_group, 'small', tables.Float32Atom(), shape=dset_small_shape, filters=compression_filter
        )
        print(
            'Writing resized tseries (TZHW={}) to "/imaging/small" (TCZHW={})'.format(
                tseries_small.shape, dset_small.shape
            )
        )
        assign_tseries_to_dset(dset_small, tseries_small)
        dset_small.attrs['dimensions'] = 'TCZHW'

        if stim_masks is not None and stim_used_at_each_timestep is not None:
            # Not creating chunked arrawys for stim data as, in this format, it really does not require much space at all
            stim_group = tyh5.create_group('/', 'stim')
            print(
                'Writing stim_masks (NZHW={}) and stim_used_at_each_timestep (T={}) to "/stim"'.format(
                    stim_masks.shape, stim_used_at_each_timestep.shape[0]
                )
            )
            mask_group = tyh5.create_group(stim_group, 'masks')
            tyh5.create_array(mask_group, 'raw', stim_masks, atom=tables.UInt8Atom())
            stim_masks_small = resize_images(stim_masks, small_size, small_size, interpolation=cv2.INTER_NEAREST)
            tyh5.create_array(mask_group, 'small', stim_masks_small, atom=tables.UInt8Atom())
            tyh5.create_array(
                stim_group, 'stim_used_at_each_timestep', stim_used_at_each_timestep, atom=tables.Int64Atom()
            )
        else:
            print('"/imaging/stim" not populated as stim data was not provided')

# https://matthew-brett.github.io/teaching/mutual_information.html
def mutual_information(A, B):
    """ Mutual information for joint histogram.
    """
    hgram, x_edges, y_edges = np.histogram2d(A.ravel(), B.ravel(), bins=20) 
    # Convert bins counts to probability values
    pxy = hgram / float(np.sum(hgram))
    px = np.sum(pxy, axis=1) # marginal for x over y
    py = np.sum(pxy, axis=0) # marginal for y over x
    px_py = px[:, None] * py[None, :] # Broadcast to multiply marginals
    # Now we can do the calculation using the pxy, px_py 2D arrays
    nzs = pxy > 0 # Only non-zero pxy values contribute to the sum
    return np.sum(pxy[nzs] * np.log(pxy[nzs] / px_py[nzs]))

def transparent_cmap(cmap, N=255, max_alpha=0.8):
    "Copy colormap and set alpha values"

    mycmap = cmap
    mycmap._init()
    mycmap._lut[:,-1] = np.linspace(0, max_alpha, N+4)
    return mycmap