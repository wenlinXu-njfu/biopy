#!/usr/bin/env python
"""
File: show_info.py
Description: Show other information.
CreateDate: 2023/8/4
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import List
from click import echo


class Displayer:
    def __init__(self,
                 program: str,
                 authors: List[str] = ['Wenlin Xu'],
                 contacts: List[str] = ['wenlinxu.njfu@outlook.com'],
                 version: str = '1.1.0'):
        self.program = program
        self.authors = ', '.join(authors)
        self.contacts = ', '.join(contacts)
        self.version = version

    def version_info(self, ctx, param, value):
        if not value or ctx.resilient_parsing:
            return
        echo(f'Program: {self.program}\n'
             f'Author:  {self.authors}\n'
             f'Contact: {self.contacts}\n'
             f'Version: {self.version}',
             err=True)
        ctx.exit()
