from os import system, listdir, remove
from os.path import abspath, isfile
from shutil import which

path = '/'.join(abspath(__file__).split('/')[:-1])
system(command="pip install fire pybioinformatic==0.0.8b  requests scipy venn")
system(command=f"sed -i 's/\r//' {path}/bin/*")
python_path = which('python')
system(command=f"sed -i 's@#!/usr/bin/env python@#!{python_path}@' {path}/bin/*")
system(command=f"chmod 755 {path}/bin/*")
system(command=f"python {path}/bin/batch_rename.py -d {path}/bin -old .py")
for file in listdir(path):
    if isfile(path + f'/{file}'):
        remove(path + f'/{file}')
