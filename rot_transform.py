from scipy.signal import fftconvolve

def apply_rot_transform(transf_data, filters):
    """ Apply rotated filters to ``data``. Contrarily to :func:`apply_transform`,
    the number of filters should match the number of scales of ``transf_data``
    (in the second axis)
    
    Parameters
    ----------
    transf_data: np.ndarray
        Transformed data. Should be in regular format (first/0th axis is image
        index), with second/1st axis indexing scales per object.
    filters: np.ndarray
        Set of (likely rotated!) filters.
    """
    rot_data = np.array([[fftconvolve(scale, filt, mode='same') for scale,filt 
                          in zip(im, filters)] for im in transf_data])
    return rot_data


