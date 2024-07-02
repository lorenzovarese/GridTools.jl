# GridTools

[![Build Status](https://github.com/jeffzwe/GridTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jeffzwe/GridTools.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Static Badge](https://img.shields.io/badge/docs-stable-blue.svg)](https://jeffzwe.github.io/GridTools.jl/dev)

## Installation

### Setup python virtual environment


### Development installation

```bash
export GRIDTOOLS_JL_PATH="..."
export GT4PY_PATH="..."
# create python virtual environemnt
#  make sure to use a python version that is compatible with GT4Py
python -m venv venv
# activate virtual env
#  this command has be run everytime GridTools.jl is used
source venv/bin/activate
# clone gt4py
git clone --branch fix_python_interp_path_in_cmake git@github.com:tehrengruber/gt4py.git
#git clone git@github.com:GridTools/gt4py.git $GT4PY_PATH
pip install -r $GT4PY_PATH/requirements-dev.txt
pip install -e $GT4PY_PATH
# 
```

## Troubleshooting

__undefined symbol: PyObject_Vectorcall__

Make sure to run everything in the same environment that you have build `PyCall` with. A common reason is you have built PyCall in a virtual environement and then didn't load it when executing stencils.