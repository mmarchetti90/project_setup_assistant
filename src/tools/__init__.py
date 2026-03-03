#!/usr/bin/env python3

from os import listdir
from os.path import abspath, dirname

init_dir = dirname(abspath(__file__))

modules = [file.replace('.py', '') for file in listdir(init_dir) if file.endswith('.py') and file != '__init__.py']

__all__ = modules
