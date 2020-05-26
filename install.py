#/////////////////////////////////////////////////////////////////////////////#
#                                                                             #
#   Copyright (C) 2020 by David B. Blumenthal                                 #
#                                                                             #
#   This file is part of GenEpiSeeker.                                        #
#                                                                             #
#   GenEpiSeeker is free software: you can redistribute it and/or modify it   #
#   under the terms of the GNU General Public License as published by         #
#   the Free Software Foundation, either version 3 of the License, or         #
#   (at your option) any later version.                                       #
#                                                                             #
#   GenEpiSeeker is distributed in the hope that it will be useful,           #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              #
#   GNU General Public License for more details.                              #
#                                                                             #
#   You should have received a copy of the GNU General Public License         #
#   along with GenEpiSeeker. If not, see <http://www.gnu.org/licenses/>.      #
#                                                                             #
#/////////////////////////////////////////////////////////////////////////////#

##
# @file install.py
# @brief Installs GenEpiSeeker and its dependencies.
#
# @details 
# Usage: 
# ```sh
# $ python install.py [--help] [-h] [--target docs|unit|instance|variance_model|penetrance_model|regression_model|bayesian_model|compare_models] [--clean] [--debug]
# ```
#
# For more information, execute `$ python install.py --help`.

'''Installs GenEpiSeeker and its dependencies.'''

import subprocess
import argparse
import os.path
import platform

def build_targets(args):  
    
    if args.clean:
        print("-- Cleaning build directory.")
        commands = "rm -rf build"
        subprocess.call(commands, shell=True)
    
    print("-- Creating build directory.")
    commands = "mkdir -p build"
    subprocess.call(commands, shell=True)
    
    if (not os.path.isfile("build/Makefile")):
        print("-- Running CMake.")
        commands = "cd build; rm -rf *; cmake .. -DCMAKE_BUILD_TYPE="
        if args.debug:
            commands = commands + "Debug"
        else:
            commands = commands + "Release"
        if platform.system() == "Darwin":
            commands = commands + " -DOMP_ROOT=" + subprocess.check_output("brew --prefix", shell=True).decode("utf-8")
        subprocess.call(commands, shell=True)
        
    print("-- Building target {}.".format(args.target))
    commands = "cd build; make {}".format(args.target)
    subprocess.call(commands, shell=True)
    
    print("-- Creating output directories.")
    commands = "mkdir -p test/unit/init; mkdir -p test/model/res"
    subprocess.call(commands, shell=True)

def build_external_libraries():
    if os.path.isfile("ext/boost_1_71_0/.INSTALLED"):
        print("-- Boost libraries already built.")
    else:
        print("-- Building Boost libraries.")
        commands = "cd ext/boost_1_71_0; ./bootstrap.sh; ./b2"
        subprocess.call(commands, shell=True)
        f = open("ext/boost_1_71_0/.INSTALLED", "w")
        f.close()

print("\n**************************************************")
print("                   GenEpiSeeker                   ")
print("                Installation Script               ")
print("**************************************************")

parser = argparse.ArgumentParser(description="Compiles GenEpiSeeker and its dependencies.", epilog="If called without arguments, only the dependencies are installed.")
targets = ["docs", "unit", "instance", "variance_model", "penetrance_model", "regression_model", "bayesian_model", "compare_models"]
parser.add_argument("--target", help="build selected target", metavar="docs|unit|instance|variance_model|penetrance_model|regression_model|bayesian_model|compare_models", choices=targets)
parser.add_argument("--debug", help="build in debug mode", action="store_true")
parser.add_argument("--clean", help="clean build directory and update makefile before build", action="store_true")
args = parser.parse_args()
build_external_libraries()
if args.target:
    build_targets(args)
print("**************************************************\n")
