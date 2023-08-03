#!/usr/bin/env python
"""
File: id_index_seq.py
Description: Extract sequences from FASTA file based on the id provided.
Date: 2022/3/25
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib.fasta import Fasta


def main(fasta_files,
         parse_seqids: bool,
         id_file,
         match: bool,
         log,
         out_file):
    content = []
    id_list = list(set(line.strip() for line in id_file))
    id_list_index = list(range(len(id_list)))  # Create a subscript for each ID
    for fasta_file in fasta_files:
        fasta_file = str(fasta_file)
        file_name = fasta_file.split('/')[-1]
        for seq_obj in Fasta(fasta_file).parse(parse_seqids):
            if match and seq_obj.id in id_list:
                del id_list[id_list.index(seq_obj.id)]
                if out_file:
                    content.append(f">{seq_obj.id} from={file_name}\n{seq_obj.seq}\n")
                else:
                    seq_obj.id = f'{seq_obj.id} from={file_name}'
                    print(seq_obj)
            elif not match:
                for _id in id_list:
                    if seq_obj.id in _id or _id in seq_obj.id:
                        index = id_list.index(_id)
                        id_list_index[index] = True  # Mark sequence that has been found
                        if out_file:
                            content.append(f">{seq_obj.id} from={file_name}\n{seq_obj.seq}\n")
                        else:
                            seq_obj.id = f'{seq_obj.id} from={file_name}'
                            print(seq_obj)
    if content:
        with open(out_file, 'w') as o:
            o.write(''.join(content))
    # report sequence that not match
    if match and id_list:
        msg = ' not found\n'.join(id_list) + ' not found\n'
        click.echo(f'\033[33m{msg}\033[0m', err=True, file=open(log, 'a')) if log else \
            click.echo(f'\033[33m{msg}\033[0m', err=True)
    # report sequence that not contain
    if not match:
        not_found = [i for i in id_list_index if i is not True]
        not_found = [id_list[i] for i in not_found]
        if not_found:
            msg = ' not found\n'.join(not_found) + ' not found\n'
            click.echo(f'\033[33m{msg}\033[0m', err=True, file=open(log, 'a')) if log else \
                click.echo(f'\033[33m{msg}\033[0m', err=True)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('fasta_files', nargs=-1)
@click.option('-P', '--parse_seqids', 'parse_seqids', is_flag=True, flag_value=True, help='Parse sequence id.')
@click.option('-I', '--id_file', 'id_file', type=click.File('r'), help='Input id TEXT file (one id per line).')
@click.option('--match/--contain', default=True, show_default=True,
              help='Whether the id supplied should exactly match the ID of the sequence.')
@click.option('-l', '--log_file', 'log_file',
              help='Output log file, if not specified, the log will print to terminal as stderr.')
@click.option('-o', '--output_file', 'outfile',
              help='Output file, if not specified, print results to terminal as stdout.')
def run(fasta_files, parse_seqids, id_file, match, log_file, outfile):
    """Extract sequences from FASTA file based on the id provided."""
    main(fasta_files, parse_seqids, id_file, match, log_file, outfile)


if __name__ == '__main__':
    run()
