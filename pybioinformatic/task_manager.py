#!/usr/bin/env python
"""
File: task_manager.py
Description: Instance a TaskManager.
CreateDate: 2023/9/8
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union, Iterable, Callable
from io import TextIOWrapper
from subprocess import run
from datetime import datetime
from getpass import getuser
from socket import gethostname
from multiprocessing import Pool
from click import echo


class TaskManager:
    def __init__(self,
                 commands: Iterable[str] = None,
                 num_processing: int = None,
                 params: Iterable[tuple] = None,
                 log_file: Union[str, TextIOWrapper] = None):
        try:
            self.task = list(commands)
        except TypeError:
            self.task = []
        self.params = params
        self.num_processing = num_processing
        if isinstance(log_file, str):
            self.loger = open(log_file, 'w')
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
        """
        NOTICE: This method cannot run normally when it is called by multiprocessing.
                Because different subprocess will share same TaskManager.loger attribution.
        """
        echo(f'\033[33m[{getuser()}@{gethostname()}: '
             f'{datetime.now().replace(microsecond=0)}]\n$ '
             f'\033[0m\033[36m{cmd}\033[0m', self.loger, err=True)
        run(cmd, shell=True, executable="/bin/bash")

    def serial_run_cmd(self):
        for cmd in self.task:
            echo(f'\033[33m[{getuser()}@{gethostname()}: {datetime.now().replace(microsecond=0)}]\n'
                 f'$ \033[0m\033[36m{cmd}\033[0m', self.loger, err=True)
            run(cmd, shell=True, executable="/bin/bash")

    def parallel_run_cmd(self):
        if not self.task:
            echo('\033[31mError: TaskManager has no task.\033[0m', err=True)
            exit()
        pool = Pool(self.num_processing)
        for cmd in self.task:
            pool.apply_async(self.echo_and_exec_cmd, args=(cmd,))
        pool.close()
        pool.join()

    def parallel_run_func(self, func: Callable, call_back_func: Callable = None):
        results = []
        with Pool(self.num_processing) as pool:
            for param in self.params:
                ret = pool.apply_async(func, args=param, callback=call_back_func)
                results.append(ret)
            pool.close()
            pool.join()
        return results


if __name__ == '__main__':
    tkm = TaskManager()
    cmds = [f'echo Hello World x {i}' for i in range(10)]
    tkm.add_task(cmds)
    tkm.serial_run_cmd()
