"""
File: timer.py
Description: Instantiate a timer class object
Date: 2021/12/31
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from datetime import datetime


class Timer:
    def __init__(self, message: str = None):
        self.message = message

    def __call__(self, function):
        def wrapper(*args, **kwargs):
            start_time = datetime.now().replace(microsecond=0)
            if self.message:
                click.echo(f"[{datetime.now().replace(microsecond=0)}] {self.message}", err=True)
            else:
                msg = function.__name__.replace('_', ' ')
                click.echo(f"[{datetime.now().replace(microsecond=0)}] Start run {msg}.", err=True)
            end_time = datetime.now().replace(microsecond=0)
            click.echo(f'[{datetime.now().replace(microsecond=0)}] Finish in {end_time - start_time}.', err=True)
            return function(*args, **kwargs)
        return wrapper
