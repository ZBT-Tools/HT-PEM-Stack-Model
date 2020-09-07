import os
import sys
from cx_Freeze import setup, Executable

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
if sys.platform == "win32":
    base = "Win32GUI"


PYTHON_INSTALL_DIR = os.path.dirname(os.path.dirname(os.__file__))
os.environ['TCL_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tcl8.6')
os.environ['TK_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tk8.6')
icon_file = r'C:\Users\lukas\PycharmProjects\PEMFCModel\pemfc\logo-zbt.ico'
build_exe_options = {'packages': ['numpy'],
                     'includes': ['matplotlib.backends.backend_tkagg'],
                     'include_files': [(os.path.join(PYTHON_INSTALL_DIR, 'DLLs', 'tcl86t.dll'),
                                        os.path.join('lib', 'tcl86t.dll')),
                                       (os.path.join(PYTHON_INSTALL_DIR, 'DLLs', 'tk86t.dll'),
                                        os.path.join('lib', 'tk86t.dll'))
                                       # add here further files which need to be included as described in 1.
                                      ]}
packages = ["os", "numpy", "scipy", "matplotlib", "tkinter"]
options = {
    'build_exe': {
        'packages': packages,
        'includes': ['matplotlib.backends.backend_tkagg'],
        'include_files': [(os.path.join(PYTHON_INSTALL_DIR, 'DLLs',
                                        'tcl86t.dll'),
                           os.path.join('lib', 'tcl86t.dll')),
                          (os.path.join(PYTHON_INSTALL_DIR, 'DLLs',
                                        'tk86t.dll'),
                           os.path.join('lib', 'tk86t.dll')),
                          icon_file]
                 }
          }
setup(name="pemfc",
      version="0.1",
      description="PEMFC Model",
      options=options,
      executables=[Executable("app.py", base=base, icon=icon_file)])
