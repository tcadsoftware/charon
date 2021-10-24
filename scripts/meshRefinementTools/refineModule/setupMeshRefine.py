
from setuptools import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize


ext = Extension(name="meshRefine",
                sources=["meshRefine.pyx", "refine_external.cpp","cellCreator.cpp","lusolve.cpp","scalarFunction.cpp","surfaceFinder.cpp","functionDefinitions.cpp","isolineFinder.cpp","distanceFunctions.cpp","transform.cpp"],
                language="c++",cmdclass={'build_ext':build_ext},options={'-O3'})

setup(ext_modules = cythonize(ext))


