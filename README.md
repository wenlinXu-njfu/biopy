# A tiny bioinformatics and visual tool.

## Dependency of python packages
**python ≥ 3.8.5<br />
click ≥ 8.1.7<br />
natsort ≥ 8.4.0<br />
matplotlib ≥ 3.5.1<br />
numpy ≥ 1.23.1<br />
openpyxl ≥ 3.0.9<br />
pandas ≥ 1.4.2<br />
requests ≥ 2.26.0<br />
scipy ≥ 1.9.0<br />
seaborn ≥ 0.11.2<br />
tqdm ≥ 4.62.3<br />
venn ≥ 0.1.3<br />**

## Dependency of other software
**blast+ for circRNA flanking sequence analyse.<br />
pfamscan for batch perform pfamscan.<br />**

## Getting started
**git clone https://github.com/wenlinXu-njfu/biopy_v1.1.0.git <br />
python biopy_v1.1.0/configure.py<br />
export PYTHONPATH=$PATH:/home/user/software/biopy_v1.1.0<br />
export PATH=$PATH:/home/user/software/biopy_v1.1.0/bin<br />**

## example
**plot circos \ <br /> -c biopy/plot_lib/circos/test_data/Ptc_chr_len.txt \ <br /> -d biopy/plot_lib/circos/test_data/gene_density.txt \ <br /> -s biopy/plot_lib/circos/test_data/stat.txt \ <br /> -l biopy/plot_lib/circos/test_data/link.txt \ <br /> -o biopy/plot_lib/circos/test_data/circos.png**
![image](plot_lib/circos/test_data/circos.png)
