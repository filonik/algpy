"""Execute the example notebooks with a fresh kernel for each notebook.

Run from an environment with algpy and the examples extra installed:
    python examples/execute.py
    python examples/execute.py --output-dir /tmp/algpy-executed
"""
import argparse
from pathlib import Path

import nbformat
from nbclient import NotebookClient


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("notebooks", nargs="*", type=Path,
                        help="Notebook paths; defaults to all notebooks beside this script")
    parser.add_argument("--kernel", default="python3", help="Installed Jupyter kernel name")
    parser.add_argument("--output-dir", type=Path,
                        help="Optional directory for executed copies; source notebooks stay clean")
    args = parser.parse_args()
    directory = Path(__file__).resolve().parent
    paths = args.notebooks or sorted(directory.glob("*.ipynb"))
    for path in paths:
        path = path.resolve()
        notebook = nbformat.read(path, as_version=4)
        nbformat.validate(notebook)
        print(f"Executing {path.name} ...", flush=True)
        NotebookClient(notebook, timeout=180, kernel_name=args.kernel,
                       resources={"metadata": {"path": str(path.parent)}}).execute()
        if args.output_dir:
            args.output_dir.mkdir(parents=True, exist_ok=True)
            nbformat.write(notebook, args.output_dir/path.name)
        print(f"Passed {path.name}", flush=True)


if __name__ == "__main__":
    main()
