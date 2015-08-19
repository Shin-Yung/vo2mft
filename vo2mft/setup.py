from argparse import ArgumentParser
import os

def _main():
    parser = ArgumentParser(description="Set up config")
    parser.add_argument('--base_path', type=str, default=None,
            help="vo2mft base path (above vo2mft and vo2solve directories)")
    args = parser.parse_args()

    base_path = args.base_path
    if args.base_path == None:
        print("base_path not specified; assuming it is one directory up from pwd.")
        pwd = os.getcwd()
        dir_up = os.path.join(pwd, '..')
        base_path = os.path.abspath(dir_up)

    config = {'base_path': base_path}

    config_path = os.path.join(base_path, 'vo2mft', 'config')

    with open(config_path, 'w') as fp:
        for k, v in config.items():
            fp.write("{} {}\n".format(k, v))

    print("config file written.")

if __name__ == "__main__":
    _main()
