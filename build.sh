#!/bin/bash
# This script make a source distribution for the pip installable package in the current folder
echo "Build package for exoplot"
python --version
# Explicit give format otherwise a zip is created (Windows?)
python setup.py sdist --formats=gztar

# Next command will not close the window, can be handy if something goes wrong
exec $SHELL