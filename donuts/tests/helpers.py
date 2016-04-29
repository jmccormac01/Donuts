import os
from astropy.utils.data import download_file

ROOT_URL = 'https://s3-eu-west-1.amazonaws.com/donuts-data'

def get_test_filename(filename, cache=True, timeout=30):
    '''Get a test filename from the package directory
    '''
    url = os.path.join(ROOT_URL, filename)
    return download_file(url, cache=cache, timeout=timeout)
