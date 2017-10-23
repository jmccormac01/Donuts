import os
from astropy.utils.data import download_file

ROOT_URL = 'http://www.jamesjmccormac.com/donuts/data'

def get_test_filename(filename, cache=True, timeout=30):
    '''Get a test filename from the package directory
    '''
    # this stopped working when we went to testing files from the web
    #url = os.path.join(ROOT_URL, filename)
    url = "{}/{}".format(ROOT_URL, filename)
    return download_file(url, cache=cache, timeout=timeout)
