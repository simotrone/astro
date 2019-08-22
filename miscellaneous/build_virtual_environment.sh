#!/bin/bash

# Prepare local environment:
python3 -m venv venv3

# Activate python3 environment:
source venv3/bin/activate

# Update pip and install modules:
pip install --upgrade pip
pip install astropy
# Could need tkinter on box (ex: dnf install python3-tkinter)
pip install matplotlib
