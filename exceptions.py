class Error(Exception):
    '''Base class for other exceptions'''
    pass

class PropertyNotFoundError(Exception):
    '''Raised when a property is not found on the output'''
    pass