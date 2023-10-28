#!/usr/bin/env python
"""
File: decompressing_file.py
Description: Decompressing file.
CreateDate: 2022/4/22
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from gzip import GzipFile


def ungz(gz_file):
    """Decompressing gz file."""
    with open(gz_file.replace('.gz', ''), 'ab') as out:
        for line in GzipFile(gz_file):
            out.write(line)
    return gz_file.replace('.gz', '')
