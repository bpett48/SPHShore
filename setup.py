"""Build the existing C++17 solver as a native CPython extension."""
import sys
from setuptools import setup, Extension

if sys.platform == "win32":
    flags = ["/std:c++17", "/O2", "/EHsc"]
else:
    flags = ["-std=c++17", "-O3"]

setup(ext_modules=[Extension(
    "sphshore._core",
    sources=["cpp/bindings.cpp", "cpp/wave_tank_solver.cpp"],
    language="c++",
    extra_compile_args=flags,
)])
