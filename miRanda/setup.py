from distutils.core import setup, Extension

module1 = Extension('pymiranda',
  sources = [
    'src/pymiranda.c',
    'src/scan.c',
    'src/swat.c',
    'src/thermo.c',
    'src/utils.c',
    'src/output.c',
    'src/ExpString.c',
  ],
  include_dirs = [
    '/usr/local/include/ViennaRNA'
  ],
  library_dirs = [
    '/usr/local/lib'
  ],
  libraries = [
    'RNA',
    'm'
  ]
)

setup (name = 'pymiranda',
        version = '1.0',
        description = 'This is a demo package',
        ext_modules = [module1])
