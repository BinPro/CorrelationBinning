
class LevelError(Exception):
    """Indicating erroneous as a taxonomic level """
    def __init__(self,level):
        self.level = level
    def __str__(self):
        return repr(self.level)
