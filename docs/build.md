# Build

Requires [Nix](https://nixos.org/download/). All dependencies are pinned and fetched automatically.

```bash
nix build
```

Outputs:

```
result/lib/libsemel.so.0.1.0   # shared library
result/lib/libsemel.a           # static library
result/include/semel/*.h        # headers
result/etc/semel/zlog.conf      # logging config
```

## Development shell

```bash
nix develop
```

This provides all C dependencies plus Python, numpy, and matplotlib.

To build manually inside the dev shell:

```bash
mkdir build && cd build
cmake -DB_PROJECT_NAME=semel \
      -DB_VERSION_MAJOR=0 \
      -DB_VERSION_MINOR=1 \
      -DB_VERSION_PATCH=0 \
      -DB_MACHINE=$(uname -m) ..
make
```

If qhull is installed outside the nix environment, pass its prefix:

```bash
cmake ... -DQHULL_PREFIX=/path/to/qhull ..
```

## Python setup

Inside `nix develop`, create a venv and install the bindings:

```bash
python -m venv .venv --system-site-packages
source .venv/bin/activate

export SEMEL_LIB_PATH=./result/lib
export LD_LIBRARY_PATH=./result/lib:$LD_LIBRARY_PATH
export SEMEL_ZLOG_CONF=./result/etc/semel/zlog.conf

pip install ./wrap/python
```

`--system-site-packages` gives the venv access to numpy, matplotlib, etc. provided by the Nix dev shell.

## Environment variables

| Variable | Description |
| --- | --- |
| `SEMEL_LIB_PATH` | Directory containing `libsemel.so` (for Python bindings) |
| `SEMEL_ZLOG_CONF` | Path to zlog.conf (overrides compile-time default) |

## Dependencies

All managed by the Nix flake:

| Library | Purpose |
| --- | --- |
| [qhull](http://www.qhull.org/) | Delaunay triangulation |
| [Judy](http://judy.sourceforge.net/) | Sparse dynamic arrays (simplices, filtration, cochains) |
| [FLINT](http://www.flintlib.org/) | Finite field and integer arithmetic for cohomology |
| [GMP](https://gmplib.org/) | Multi-precision arithmetic (FLINT dependency) |
| [MPFR](https://www.mpfr.org/) | Multi-precision floats (FLINT dependency) |
| [GSL](https://www.gnu.org/software/gsl/) | Linear algebra (circumsphere, Hodge Laplacian) |
| [zlog](https://github.com/HardySimpson/zlog) | Logging |

Logging is configured in `config/zlog.conf` (logs go to `/var/log/semel/`).
