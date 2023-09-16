#!/usr/bin/env python
"""
File: task_manager.py
Description: Instance a TaskManager.
Date: 2023/9/8
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union, Iterable
from io import TextIOWrapper
from os import system
from datetime import datetime
from getpass import getuser
from socket import gethostname
from multiprocessing import Pool
from click import echo


class TaskManager:
    def __init__(self,
                 commands: Iterable[str] = None,
                 processing_num: int = None,
                 log_file: Union[str, TextIOWrapper] = None):
        try:
            self.task = list(commands)
        except TypeError:
            self.task = []
        self.processing_num = processing_num
        if isinstance(log_file, str):
            self.loger = open(log_file, 'a')
        else:
            self.loger = log_file

    def add_task(self, command: Union[str, Iterable[str]]):
        if isinstance(command, str):
            self.task.append(command)
        else:
            command = list(command)
            self.task.extend(command)

    def del_task(self, index: Union[int, str] = None):
        if isinstance(index, int):
            del self.task[index]
        elif isinstance(index, str):
            self.task.remove(index)
        else:
            self.task.pop()

    def clear_task(self):
        self.task = []

    def echo_and_exec_cmd(self, cmd: str):
        echo(f'\033[33m[{getuser()}@{gethostname()}: '
             f'{datetime.now().replace(microsecond=0)}]\n$ '
             f'\033[0m\033[36m{cmd}\033[0m', self.loger, err=True)
        system(cmd)

    def serial_run(self):
        for cmd in self.task:
            echo(f'\033[33m[{getuser()}@{gethostname()}: {datetime.now().replace(microsecond=0)}]\n'
                 f'$ \033[0m\033[36m{cmd}\033[0m', self.loger, err=True)
            system(cmd)

    def parallel_run(self):
        if not self.task:
            echo('\033[31mError: TaskManager has no task.\033[0m', err=True)
            exit()
        pool = Pool(self.processing_num)
        for cmd in self.task:
            pool.apply_async(self.echo_and_exec_cmd, args=(cmd,))
        pool.close()
        pool.join()


if __name__ == '__main__':
    tkm = TaskManager()
    cmds = [f'echo Hello World x {i}' for i in range(10)]
    tkm.add_task(cmds)
    tkm.serial_run()
