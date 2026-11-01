import argparse

from pnatlas import __version__


def main(argv=None):
    parser = argparse.ArgumentParser(prog="pnatlas")
    parser.add_argument("--version", action="version", version=f"pnatlas {__version__}")
    parser.parse_args(argv)
    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
