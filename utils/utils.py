def pp(d, t=0):
    TAB = "    "
    result = ",\n".join([
        f"{TAB * (t + 1)}{k} : {v if not isinstance(v, dict) else pp(v, t+1)}"
        for k, v in d.items()
    ])

    return "{\n" + result + f"\n{TAB * t}}}"


def deep_get(dict, path, default=None):
    if not isinstance(path, tuple):
        path = (path, )

    for elem in path:
        if elem not in dict:
            return default
        dict = dict[elem]
    return dict
