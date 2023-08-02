from os import system, listdir, remove
from os.path import abspath, isfile

path = abspath(__file__).replace(__file__, '')
system(command=f"sed -i 's/\r//' {path}/bin/*")
system(command=f"chmod u+x {path}/bin/*")
system(command=f"python {path}/bin/batch_rename.py -d {path}/bin -old .py")
for file in listdir(path):
    if isfile(path + f'/{file}'):
        remove(path + f'/{file}')
