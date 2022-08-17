class ZeroDict(dict):
    """
    A dict that returns 0 for unknown keys.
    """
    def __getitem__(self, item):
        return self.get(item)

    def get(self, key):
        return super().get(key, 0)
