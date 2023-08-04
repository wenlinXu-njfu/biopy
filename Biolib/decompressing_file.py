#!/usr/bin/env python
"""
File: decompressing_file.py
Description: Decompressing file
Date: 2022/4/22
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union, IO
import gzip


def ungz(gz_file: Union[IO, str]):
    """Decompressing gz file."""
    with open(gz_file.replace('.gz', ''), 'ab') as out:
        for line in gzip.GzipFile(gz_file):
            out.write(line)
    return gz_file.replace('.gz', '')
