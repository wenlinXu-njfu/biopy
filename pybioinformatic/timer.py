"""
File: timer.py
Description: Instantiate a timer class object
Date: 2021/12/31
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from datetime import datetime
from click import echo


class Timer:
    def __init__(self, start_message: str = None, command_content: str = None):
        self.start_message = start_message
        self.command_content = command_content

    def __call__(self, function):
        def wrapper(*args, **kwargs):
            start_time = datetime.now().replace(microsecond=0)
            if self.start_message:
                echo(f"[{datetime.now().replace(microsecond=0)}] {self.start_message}", err=True)
            else:
                msg = function.__name__.replace('_', ' ')
                echo(f"[{datetime.now().replace(microsecond=0)}] Start run {msg}.", err=True)
            if self.command_content:
                echo(f"[{datetime.now().replace(microsecond=0)}] {self.command_content}", err=True)
            value = function(*args, **kwargs)
            end_time = datetime.now().replace(microsecond=0)
            echo(f'[{datetime.now().replace(microsecond=0)}] Finish in {end_time - start_time}.', err=True)
            return value
        return wrapper