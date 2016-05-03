class Image(object):
    '''Encapsulate the transformations applied to images
    '''
    def __init__(self, data, header):
        self.raw_image = data
        self.header = header
