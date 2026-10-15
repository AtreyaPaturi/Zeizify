import sys

from .cli import main

if __name__ == "__main__":  # guard: analysis worker processes re-import this module on Windows
    sys.exit(main())
