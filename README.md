# GridTools

[![Build Status](https://github.com/jeffzwe/GridTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jeffzwe/GridTools.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Static Badge](https://img.shields.io/badge/docs-stable-blue.svg)](https://jeffzwe.github.io/GridTools.jl/dev)

## Installation

### Development Installation

Follow these steps for setting up the development environment for `GridTools.jl`.

1. **Set Environment Variables**
   Define the paths for `GridTools.jl` and `GT4Py`.

   ```bash
   export GRIDTOOLS_JL_PATH="/path/to/GridTools.jl"
   export GT4PY_PATH="/path/to/gt4py"
   ```

2. **Create a Python Virtual Environment**
   Ensure you're using a Python version that is compatible with `GT4Py`.

   ```bash
   python -m venv .venv
   ```

3. **Activate the Virtual Environment**
   This command needs to be run every time you use `GridTools.jl`.

   ```bash
   source .venv/bin/activate
   ```

4. **Clone `gt4py` Repository**
   Clone the `gt4py` repository and switch to the required branch.

   ```bash
   git clone --branch fix_python_interp_path_in_cmake git@github.com:tehrengruber/gt4py.git $GT4PY_PATH
   # git clone git@github.com:GridTools/gt4py.git $GT4PY_PATH
   ```

5. **Install Dependencies**
   Install the development requirements and `gt4py` in editable mode.

   ```bash
   pip install -r $GT4PY_PATH/requirements-dev.txt
   pip install -e $GT4PY_PATH
   ```

## Troubleshooting

__undefined symbol: PyObject_Vectorcall__

Make sure to run everything in the same environment that you have build `PyCall` with. A common reason is you have built PyCall in a virtual environement and then didn't load it when executing stencils.

### Python Versions and PyCall

If you use a Python version greater than 3.10.14, you might encounter some errors. Follow these steps for the correct execution of `GridTools.jl`.

1. Use/Activate the correct python version, manually or with pyenv:
```bash
# Example with pyenv
pyenv global 3.10.14
pyenv rehash
```

2. Create a virtual enviroment
```bash
python -m venv .venv
```

3. Activate the virtual environment
```bash
source .venv/bin/activate
```

4. Set the PYTHON environment variable for Julia
```bash
export PYTHON="/path/to/.venv/bin/python"
```

5. Rebuild PyCall with the correct Python
```bash
julia -e 'using Pkg; ENV["PYTHON"] = "/path/to/.venv/bin/python"; Pkg.build("PyCall")'
```

6. Verify the correctness of the configuration running GridTools.jl
```bash
# Run the Julia script
julia --color=yes --project=/path/to/.julia/environments/v1.10 /path/to/GridTools.jl/src/GridTools.jl
```
