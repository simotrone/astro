#!/bin/bash

# Prepare local environment:
python3 -m venv venv3

# Activate python3 environment:
source venv3/bin/activate

# Update pip and install modules:
pip install --upgrade pip
pip install astropy
