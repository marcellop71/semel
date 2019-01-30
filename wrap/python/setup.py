from setuptools import setup
from setuptools import find_packages

setup(name='semel',
      version='0.1.0',
      description='Simplicial methods library',
      author='Marcello Paris',
      author_email='marcello.paris@gmail.com',
      license='MIT',
      install_requires=['numpy'],
      extras_require={
          'viz': ['matplotlib'],
          'hodge': ['torch'],
          'transport': ['pot'],
      },
      packages=find_packages())
