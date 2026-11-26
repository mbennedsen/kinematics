# Native SHA-256 core

`sha256Compress.m`, `sha256Finalize.m` and `sha256Initialize.m` are
vendored from [DataHash](https://github.com/stefanopradella/DataHash),
commit `21d6e047430e0cab7602dba23759ed20e1f2312b` (retrieved 26 November 2026).
The upstream files are under `src/sha256/`.

These files use the upstream algorithms, with attribution headers and a
licence reference to the bundled `DataHash_LICENSE.txt` (BSD 3-Clause terms).

The package's `file_sha256.m` supplies its own bounded-memory binary file
reader and hexadecimal output wrapper. All checksum calculations run in MATLAB.
