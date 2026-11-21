"""Lets `python -m refmasker ...` work the same as the `refmasker` command."""

from refmasker.cli import main

if __name__ == "__main__":
    # main() returns the exit code; SystemExit hands it to the shell.
    raise SystemExit(main())
