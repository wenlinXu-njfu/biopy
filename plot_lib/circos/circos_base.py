"""
File: circos_base.py
Description: Get theta angle and width of each chromosome.
CreateDate: 2022/10/10
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import numpy as np
from typing import Tuple


def get_chr_theta_width(chr_len_file: str) -> Tuple[dict, dict, dict]:
    """
    Get theta angle and width of each chromosome bar.
    :param chr_len_file: Chromosome length file. (Chr_num\tlength\n)
    :return: tuple(raw_chr_len_dict: dict, chr_theta_dict: dict, chr_width_dict: dict)
             raw_chr_len_dict -> raw chromosome length dict. {Chr_num: length, ...}
             chr_theta_dict -> chromosome theta dict. {Chr_num: theta, ...}
             chr_width_dict -> chromosome width dict. {Chr_num: width, ...}
    """
    # Step 1: Read in chromosome length file.
    raw_chr_len_dict = {line.split('\t')[0]: int(line.strip().split('\t')[1])
                        for line in open(chr_len_file)}
    # Step 2: Calculate length percentage of each chromosome.
    chr_width_dict = {}
    each_space_len = min(raw_chr_len_dict.values()) / 6
    total_space_len = each_space_len * len(raw_chr_len_dict.keys())
    for chr_num, chr_len in raw_chr_len_dict.items():
        chr_width_dict[chr_num] = 2 * np.pi * chr_len / (sum(raw_chr_len_dict.values()) + total_space_len)
    # Step 3: Calculate theta angle of each chromosome.
    chr_theta_dict = {}
    chr_theta = 0
    space_angle = 2 * np.pi * each_space_len / (sum(raw_chr_len_dict.values()) + total_space_len)
    for chr_num, chr_len in chr_width_dict.items():
        if chr_theta == 0:
            chr_theta_dict[chr_num] = chr_len / 2
        else:
            chr_theta_dict[chr_num] = chr_theta + chr_len / 2
        chr_theta += chr_len + space_angle
    return raw_chr_len_dict, chr_theta_dict, chr_width_dict
