import cv2
import numpy as np
import tables
from tqdm import tqdm


def assign_tseries_to_dset(dset, tseries, channel=0, batch_size=1024):
    "dset is 5D (T, C, Z, H, W), tseries is 4D (T, Z, H, W)"
    if len(tseries.shape) == 4:
        T, Z, H, W = tseries.shape
        for z in range(Z):
            print(f"saving z-plane {z}")
            # Using batching to write out tseries (these blobs are massive ~75GB)
            for b in tqdm(range(0, T, batch_size)):
                dset[b:b + batch_size, channel, z] = tseries[b:b + batch_size, z]
    else:
        raise NotImplementedError


def resize_images(images, target_width, target_height, interpolation=cv2.INTER_LINEAR):
    assert len(images.shape) == 4
    T, Z, original_height, original_width = images.shape
    print('Original: {}'.format(images.shape))
    fx = target_width / original_width
    fy = target_height / original_height

    example_resized_image = cv2.resize((images[0, 0]), None, fx=fx, fy=fy, interpolation=interpolation)
    assert example_resized_image.shape == (target_height, target_width), 'Only allowing exact resizes'

    resized_images = np.zeros((images.shape[0], images.shape[1], target_height, target_width), dtype=np.float32)
    for t, volume in enumerate(images):
        for z, image in enumerate(volume):
            resized_images[t, z] = cv2.resize(image, None, fx=fx, fy=fy, interpolation=interpolation)
    return resized_images


def write_tseries_to_tyh5(tseries, output_path, num_channels=1, compression_level=3, small_size=256):
    """
    NOTE: not a huge deal if doesn't fit in memory (can just read from tiff_files)

    Parameters
    ----------
    tseries

    output_path

    num_channels

    compression_level

    Notes
    -----
    This function is less useful if full tseries doesn't fit in RAM (in that case we'll
    need to load from TIFF files directly - not the end of the world)
    """
    print(tseries.shape)
    print(tseries.dtype)
    print('Output path: {}'.format(output_path))
    with tables.open_file(output_path, 'w') as tyh5:
        print(f'opened {output_path}')

        # This doesn't seem to work
        # if not 'version' in tyh5.root._v_attrs:
        #     tyh5.root._v_attrs['version'] = version

        compression_filter = tables.filters.Filters(complevel=compression_level, complib='blosc:zstd')

        # First just assume shape that comes out of Lensman.loadTseries
        # TODO(allan.raventos): add support for channels if needed
        assert len(tseries.shape) == 4
        H, W, Z, T = tseries.shape
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
        print('tseries shape before adding {}'.format(tseries.shape))
        assign_tseries_to_dset(dset_raw, tseries)
        dset_raw.attrs['dimensions'] = 'TCZHW'

        # Need to generate 256. TODO(allan.raventos): generate 512 if raw is 1024
        # TB comment: perhaps convenient if small is always 256, regardless of raw
        dset_small_shape = (T, num_channels, Z, small_size, small_size)
        print('about to resize {}'.format(dset_small_shape))
        tseries_small = resize_images(tseries, small_size, small_size, interpolation=cv2.INTER_LINEAR)
        dset_small = tyh5.create_carray(
            imaging_group, 'small', tables.Float32Atom(), shape=dset_small_shape, filters=compression_filter
        )
        assign_tseries_to_dset(dset_small, tseries_small)
        dset_small.attrs['dimensions'] = 'TCZHW'

        # TODO(allan.raventos): write out stim mask
        # Tensor holding masks (num_masks, H, W)
        # 1D array with entry i corresponding to mask used at timestep i (could use 0 for no mask)
