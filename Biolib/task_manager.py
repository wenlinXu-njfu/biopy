#!/usr/bin/env python
"""
File: task_manager.py
Description: Instance a TaskManager.
Date: 2023/9/8
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Iterable
from io import TextIOWrapper
from os import system
from datetime import datetime
from getpass import getuser
from socket import gethostname
from multiprocessing import Pool
from click import echo


class TaskManager:
    def __init__(self, commands: Iterable[str], processing_num: int, log_file: TextIOWrapper = None):
        self.cmd_prompt = f'[{getuser()}@{gethostname()}: {datetime.now().replace(microsecond=0)}]\n$ '
        self.task = commands
        self.processing_num = processing_num
        self.loger = log_file

    def serial_run(self):
        for cmd in self.task:
            echo(f'\033[33m{self.cmd_prompt}\033[0m\033[36m{cmd}\033[0m', self.loger, err=True)
            system(cmd)

    def parallel_run(self):
        with Pool(self.processing_num) as pool:
            for cmd in self.task:
                pool.apply_async(system, (cmd,),
                                 echo(f'\033[33m{self.cmd_prompt}\033[0m\033[36m{cmd}\033[0m', self.loger, err=True))
            pool.close()
            pool.join()
