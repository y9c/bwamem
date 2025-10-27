#!/usr/bin/env python3
"""
Minimal setup.py for C extensions with uv.
All package metadata is defined in pyproject.toml.
"""

import os
import subprocess
import sys
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext

class build_ext(_build_ext):
    def build_extensions(self):
        # Build the static library first
        self.build_static_library()
        
        # Remove '-Wl,--exclude-libs,ALL' so symbols from static libs are exported
        comp = getattr(self, 'compiler', None)
        if comp and hasattr(comp, 'linker_so'):
            comp.linker_so = [
                arg for arg in comp.linker_so
                if not (isinstance(arg, str) and '--exclude-libs' in arg)
            ]
        super().build_extensions()
    
    def build_static_library(self):
        """Build the static library if it doesn't exist or is outdated."""
        static_lib_path = os.path.join('bwa', 'libbwa.a')
        
        if not os.path.exists(static_lib_path):
            print("Building static library...")
            try:
                subprocess.run(['make', 'bwa/libbwa.a'], 
                              cwd=os.getcwd(), 
                              check=True, 
                              capture_output=True, 
                              text=True)
                print("Static library built successfully")
            except subprocess.CalledProcessError as e:
                print(f"Error building static library: {e}")
                print(f"stdout: {e.stdout}")
                print(f"stderr: {e.stderr}")
                sys.exit(1)
        else:
            print("Static library already exists")

# C Extension definition
extensions = [
    Extension(
        'bwalib',
        sources=['bwamem/libbwamem.c', 'bwamem/memopts.c'],
        include_dirs=['bwa'],
        extra_compile_args=['-pedantic', '-Wall', '-std=c99', '-march=native', '-ffast-math', '-DUSE_SSE2', '-DNDEBUG'],
        extra_link_args=['-Wl,--whole-archive', os.path.join('bwa', 'libbwa.a'), '-Wl,--no-whole-archive', '-lz']
    )
]

# Minimal setup - all metadata comes from pyproject.toml
setup(
    ext_modules=extensions,
    cmdclass={'build_ext': build_ext}
)