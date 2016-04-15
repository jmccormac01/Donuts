import os
from astropy.utils.data import download_file

ROOT_URL = 'http://deneb.astro.warwick.ac.uk/phsnag/donuts'

def get_test_filename(filename, cache=True):
    '''Get a test filename from the package directory
    '''
    url = os.path.join(ROOT_URL, filename)
    return download_file(url, cache=cache)
