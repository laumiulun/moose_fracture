#!/usr/bin/env python
import sys
import os

# Locate MOOSE directory
os.chdir(os.path.abspath(os.path.dirname(__file__)))
MOOSE_DIR = os.getenv('MOOSE_DIR', os.path.join(os.getcwd(), '..', '..', 'moose'))

if not os.path.exists(MOOSE_DIR):
    raise Exception('Failed to locate MOOSE, specify the MOOSE_DIR environment variable.')

# Append MOOSE python directory
MOOSE_PYTHON_DIR = os.path.join(MOOSE_DIR, 'python')
if MOOSE_PYTHON_DIR not in sys.path:
    sys.path.append(MOOSE_PYTHON_DIR)
import MooseDocs
from MooseDocs import main
MooseDocs.MOOSE_DIR = MOOSE_DIR

if __name__ == '__main__':
    sys.exit(main.run())
