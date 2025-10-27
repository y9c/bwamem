#!/usr/bin/env python3
"""
Setup script for bwamem package.
This is a minimal setup.py that works with pyproject.toml for C extensions.
"""

import os
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext

# Define the C extension
class build_ext(_build_ext):
    def build_extensions(self):
        # Remove '-Wl,--exclude-libs,ALL' so symbols from static libs are exported
        comp = getattr(self, 'compiler', None)
        if comp and hasattr(comp, 'linker_so'):
            comp.linker_so = [
                arg for arg in comp.linker_so
                if not (isinstance(arg, str) and '--exclude-libs' in arg)
            ]
        super().build_extensions()


extensions = [
    Extension(
        'bwalib',
        sources=['bwamem/libbwamem.c', 'bwamem/memopts.c'],
        include_dirs=['bwa'],
        extra_compile_args=['-pedantic', '-Wall', '-std=c99', '-march=native', '-ffast-math', '-DUSE_SSE2', '-DNDEBUG'],
        libraries=['z'],
        extra_link_args=['-Wl,--whole-archive', os.path.join('bwa', 'libbwa.a'), '-Wl,--no-whole-archive']
    )
]

if __name__ == '__main__':
    setup(ext_modules=extensions, cmdclass={'build_ext': build_ext})