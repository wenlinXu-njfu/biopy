#!/usr/bin/env python
"""
File: arithmometer.py
Description: A arithmometer.
CreateDate: 2023/11/30
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Tuple
from fire import Fire


def SUM(*args: Tuple[float]):
    """Computes the sum of an array."""
    return sum(args)


def SUB(*args: Tuple[float]):
    """Calculate the difference of an array."""
    minuend = args[0]
    for i in args[1:]:
        minuend -= i
    return minuend


def MUL(*args: Tuple[float]):
    """Computes the product of an array."""
    factor = args[0]
    for i in args[1:]:
        factor *= i
    return factor


def DIV(*args: Tuple[float]):
    """Computes the quotient of an array."""
    dividend = args[0]
    for i in args[1:]:
        dividend /= i
    return dividend


def exp(arithmetic_exp: str):
    """Calculate arithmetic expressions."""
    return eval(arithmetic_exp)


if __name__ == '__main__':
    Fire(name=r'A arithmometer.')
