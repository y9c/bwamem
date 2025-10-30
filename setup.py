#!/usr/bin/env python3
"""
Minimal setup.py for C extensions with uv.
All package metadata is defined in pyproject.toml.
"""

import os
import subprocess
import sys
import platform
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext


class build_ext(_build_ext):
    def build_extensions(self):
        # Build the static library first
        self.build_static_library()

        # Remove '-Wl,--exclude-libs,ALL' so symbols from static libs are exported
        # This is only needed on Linux with GNU ld
        if platform.system() == "Linux":
            comp = getattr(self, "compiler", None)
            if comp and hasattr(comp, "linker_so"):
                comp.linker_so = [
                    arg
                    for arg in comp.linker_so
                    if not (isinstance(arg, str) and "--exclude-libs" in arg)
                ]
        super().build_extensions()

    def build_static_library(self):
        """Build the static library if it doesn't exist or is outdated."""
        static_lib_path = os.path.join("bwa", "libbwa.a")

        if not os.path.exists(static_lib_path):
            print("Building static library...")
            try:
                subprocess.run(
                    ["make", "bwa/libbwa.a"],
                    cwd=os.getcwd(),
                    check=True,
                    capture_output=True,
                    text=True,
                )
                print("Static library built successfully")
            except subprocess.CalledProcessError as e:
                print(f"Error building static library: {e}")
                print(f"stdout: {e.stdout}")
                print(f"stderr: {e.stderr}")
                sys.exit(1)
        else:
            print("Static library already exists")


# Determine platform-specific compilation and linking flags
def get_compile_args():
    """Get platform-specific compilation arguments."""
    args = [
        "-pedantic",
        "-Wall",
        "-std=c99",
        "-O3",
        "-fno-finite-math-only",
        "-funroll-loops",
        "-DNDEBUG",
    ]
    
    # Add architecture-specific flags only on x86-64
    machine = platform.machine().lower()
    if machine in ["x86_64", "amd64"]:
        args.append("-march=x86-64-v2")  # Modern baseline: SSE4.2
    
    return args


def get_link_args():
    """Get platform-specific linking arguments."""
    system = platform.system()
    static_lib = os.path.join("bwa", "libbwa.a")
    
    if system == "Linux":
        # Use --whole-archive on Linux to ensure all symbols are included
        return [
            "-Wl,--whole-archive",
            static_lib,
            "-Wl,--no-whole-archive",
            "-lz",
        ]
    elif system == "Darwin":
        # macOS doesn't support --whole-archive, use -force_load instead
        return [
            "-Wl,-force_load," + static_lib,
            "-lz",
        ]
    else:
        # Other platforms: just link the static library directly
        return [static_lib, "-lz"]


# C Extension definition
extensions = [
    Extension(
        "bwalib",
        sources=["bwamem/libbwamem.c", "bwamem/memopts.c"],
        include_dirs=["bwa"],
        extra_compile_args=get_compile_args(),
        extra_link_args=get_link_args(),
    )
]


# Minimal setup - all metadata comes from pyproject.toml
setup(ext_modules=extensions, cmdclass={"build_ext": build_ext})
