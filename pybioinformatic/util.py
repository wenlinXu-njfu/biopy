"""
File: util.py
Description: Util module.
CreateDate: 2024/5/12
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Callable
from natsort import natsort_key


def dict_sort_by_keys(d: dict, by: Callable = natsort_key,  reverse: bool = False) -> dict:
    return {key: d[key] for key in sorted(d, key=by, reverse=reverse)}
