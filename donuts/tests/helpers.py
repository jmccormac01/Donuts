from astropy.utils.data import get_pkg_data_filename


def get_test_filename(filename):
    '''Get a test filename from the package directory
    '''
    return get_pkg_data_filename(
        '../data/{filename}'.format(filename=filename)
    )
