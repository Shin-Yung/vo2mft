import os

def _solve_front_path():
    # TODO - avoid assumption that pwd = vo2mft/vo2mft.
    # Can use $GOPATH (in general need to consider multiple entries (separated
    # by color).
    # TODO - get from gopath/bin instead?
    config_path = "config"
    config = _load_config(config_path)
    return os.path.join(config["base_path"], "vo2solve", "vo2solve_front", "vo2solve_front")

def _load_config(config_path):
    config = {}
    with open(config_path, 'r') as fp:
        lines = fp.readlines()
        for line in lines:
            lstr = line.strip()
            first_space = lstr.index(' ')
            k = lstr[:first_space]
            v = lstr[first_space+1:]
            config[k] = v
    return config
