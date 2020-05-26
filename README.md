# GenEpiSeeker

## 1. About GenEpiSeeker

## 2. License and Citing

GenEpiSeeker is distributed under the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html).

## 3. Dependencies

Before using GenEpiSeeker, you have to install the following external dependencies:

- [CMake](https://cmake.org/) version 2.6 or higher (for compilation). Installation instructions can be found [here](https://cmake.org/install/). 
- [Doxygen](http://www.stack.nl/~dimitri/doxygen/) (for creating the documentation). Installation instructions can be found [here](https://www.stack.nl/~dimitri/doxygen/manual/install.html).
- [OpenMP](http://www.openmp.org/) compatible C++ compiler. Under Linux, OpenMP is supported by default. Under macOS, please install [libomp](https://formulae.brew.sh/formula/libomp) using [Homebrew](https://brew.sh/).  After having installed Homebrew, open a shell and execute `brew install libomp`.

The following external dependencies are distributed with GenEpiSeeker and do not have to be installed separately:

- [Boost (version 1.71.0)](https://www.boost.org/).
- [Catch (version 2.11.0)](https://github.com/catchorg/Catch2)
- [CLI11 (version 1.9.0)](https://github.com/CLIUtils/CLI11)
- [Eigen (version 3.3.7)](http://eigen.tuxfamily.org/)

## 4. Installation under Unix

After having installed the external dependencies, execute the script `install.py` for installing GenEpiSeeker:

```sh
python install.py [--help] [-h] [--target <TARGET>] [--debug] [--clean]
```

- `--help`, `-h`: Show help message and exit.  
- `--target <TARGET>`: Build selected target. Use option `--target docs` to buid the Doxygen documentation. For a list of all other available targets, print the help message.
- `--debug`: Build in debug mode. 
- `--clean`: Delete the build directoy and update the makefile before the build.

If you execute `install.py` without any arguments, only the external libraries distributed with GenEpiSeeker are installed.