#!/usr/bin/env python
"""
File: KEGG.py
Description: Batch download sequences from KEGG url.
CreateDate: 2022/6/19
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from random import choice
from time import sleep
import requests
from natsort import natsort_key
from tqdm import tqdm
import click
from pybioinformatic import TaskManager, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.3.0')


def download(url, proxy_list: list):
    header_list = [
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.99 Safari/537.36",
        "Mozilla/5.0 (Windows NT 6.3; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.95 Safari/537.36",
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_9_2) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/35.0.1916.153 Safari/537.36",
        "Mozilla/5.0 (Windows NT 6.1; WOW64; rv:30.0) Gecko/20100101 Firefox/30.0",
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_9_2) AppleWebKit/537.75.14 (KHTML, like Gecko) Version/7.0.3 Safari/537.75.14",
        "Mozilla/5.0 (compatible; MSIE 10.0; Windows NT 6.2; Win64; x64; Trident/6.0)",
        'Mozilla/5.0 (Windows; U; Windows NT 5.1; it; rv:1.8.1.11) Gecko/20071127 Firefox/2.0.0.11',
        'Opera/9.25 (Windows NT 5.1; U; en)',
        'Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1; SV1; .NET CLR 1.1.4322; .NET CLR 2.0.50727)',
        'Mozilla/5.0 (compatible; Konqueror/3.5; Linux) KHTML/3.5.5 (like Gecko) (Kubuntu)',
        'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.8.0.12) Gecko/20070731 Ubuntu/dapper-security Firefox/1.5.0.12',
        'Lynx/2.8.5rel.1 libwww-FM/2.14 SSL-MM/1.4.1 GNUTLS/1.2.9',
        "Mozilla/5.0 (X11; Linux i686) AppleWebKit/535.7 (KHTML, like Gecko) Ubuntu/11.04 Chromium/16.0.912.77 Chrome/16.0.912.77 Safari/535.7",
        "Mozilla/5.0 (X11; Ubuntu; Linux i686; rv:10.0) Gecko/20100101 Firefox/10.0"
    ]

    if proxy_list is None:
        proxy_list = [
            {'http': 'http://39.104.23.154:6379'},
            {'http': 'http://47.104.198.111:80'},
            {'http': 'http://112.124.36.55:8888'},
            {'http': 'http://47.116.210.163:8080'},
            {'http': 'http://47.119.22.92:9098'},
            {'http': 'http://39.102.214.199:8081'},
            {'http': 'http://39.102.209.128:3128'},
            {'http': 'http://8.138.82.6:8080'},
            {'http': 'http://121.43.154.123:9098'},
            {'http': 'http://39.102.214.199:3128'}
        ]

    while True:
        response = requests.get(
            url=url,
            headers={"User-Agent": choice(header_list)},
            proxies=choice(proxy_list)
        )
        if response.status_code != 403:
            break
        else:
            sleep(180)

    return response.text


def main(in_file: TextIOWrapper,
         seq_type: click.Choice(['aaseq', 'ntseq']),
         num_processing: int,
         proxy: TextIOWrapper = None,
         output_file: TextIOWrapper = None):
    id_list = list({line.strip() for line in in_file if line.strip()})
    id_list.sort(key=natsort_key)
    url_template = 'https://rest.kegg.jp/get/{dbentries}/{seq_type}'
    if proxy:
        proxy = [{line.strip().split('\t')[0]: line.strip().split('\t')[1]} for line in proxy if line.strip()]
    if len(id_list) <= 10:
        url = url_template.format(
            dbentries='+'.join(id_list),
            seq_type=seq_type
        )
        ret = download(url, proxy)
        click.echo(ret, output_file)
    else:
        params = [
            (
                url_template.format(
                    dbentries='+'.join(id_list[i:i + 10]),
                    seq_type=seq_type
                ),
                proxy
            )
            for i in range(0, len(id_list), 10)
        ]
        tkm = TaskManager(num_processing=num_processing, params=params)
        with tqdm(total=len(params)) as pbar:
            rets = tkm.parallel_run_func(func=download, call_back_func=lambda _: pbar.update())
            for ret in rets:
                click.echo(ret.get(), output_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--seq_id', 'id_file',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help=r'Input KEGG list file, one id per line. (eg. pop:112323434\npop:112323435)')
@click.option('-t', '--seq_type', 'seq_type',
              metavar='<aaseq|ntseq>', type=click.Choice(['aaseq', 'ntseq']), default='aaseq', show_default=True,
              help='Specified sequence type.')
@click.option('-n', '--num-processing', 'num_processing',
              metavar='<int>', type=int, default=10, show_default=True,
              help='Number of processing. (Download 10 sequences per processing)')
@click.option('-p', '--proxy', 'proxy',
              metavar='<file|stdin>', type=click.File('r'),
              help=r'Proxy ip. (eg. http\thttp://39.104.23.154:6379)')
@click.option('-o', '--output-file', 'output_file',
              metavar='<file|stdout>', type=click.File('w'),
              help='Output file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(id_file, seq_type, num_processing, proxy, output_file):
    """Batch download sequences from KEGG url."""
    main(id_file, seq_type, num_processing, proxy, output_file)


if __name__ == '__main__':
    run()
