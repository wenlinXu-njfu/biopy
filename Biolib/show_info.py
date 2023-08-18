#!/usr/bin/env python
"""
File: show_info.py
Description: Show other information.
Date: 2023/8/4
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from os import system
from click import echo


class Displayer:
    def __init__(self,
                 program: str,
                 author: str = 'Wenlin Xu',
                 contact: str = 'wenlinxu.njfu@outlook.com',
                 version: str = '1.1.0'):
        self.program = program
        self.author = author
        self.contact = contact
        self.version = version

    @staticmethod
    def echo_and_execute_command(command: str):
        echo(f'\033[36m{command}\033[0m', err=True)
        system(command)

    def version_info(self, ctx, param, value):
        if not value or ctx.resilient_parsing:
            return
        echo(f'Program: {self.program}\n'
             f'Author: {self.author:>10}\n'
             f'Contact: {self.contact}\n'
             f'Version: {self.version}',
             err=True)
        ctx.exit()
