def deep_get(dict, path, default=None):
    if not isinstance(path, tuple):
        path = (path, )

    for elem in path:
        if elem not in dict:
            return default
        dict = dict[elem]
    return dict
