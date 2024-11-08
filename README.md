# A tiny bioinformatics and visual tool.

## Dependence of python packages
**python ≥ 3.8<br />
fire ≥ 0.6.0<br />
pybioinformatic == 0.0.6<br />
requests ≥ 2.26.0<br />
scipy ≥ 1.9.0<br />
venn ≥ 0.1.3<br />**
## Dependency of other software
**blast+<br />
bwa<br />
cufflinks<br />
fastp<br />
featureCounts<br />
gatk<br />
hisat2<br />
hmmsearch<br />
HTseq<br />
PfamScan<br />
samtools<br />
stringtie**

## Getting started
```shell
git clone https://github.com/wenlinXu-njfu/biopy.git
python biopy/configure.py
export PYTHONPATH=$PATH:/your/path/biopy
export PATH=$PATH:/your/path/biopy/bin
```

## Issue
### ImportError: libffi.so.7: cannot open shared object file: No such file or directory
```shell
# First use the following command to verify that the file exists in that path.
ls /usr/lib/x86_64-linux-gnu/libffi.so.7

# If libffi.so.6 is present on your system but libffi.so.7 is missing, you can try creating a soft link to an existing libffi.so.6 file.
ln -s /usr/lib/x86_64-linux-gnu/libffi.so.6 /usr/lib/x86_64-linux-gnu/libffi.so.7

# You can also install libffi7 with sudo grant.
sudo apt-get install libffi7
```

## Example
### Run commands concurrently.
```shell
for i in `ls dir`;
do echo blastn -query dir/$i -db genome.fa -out "$i.blastn.xls";
done | exec_cmds -f - -n 10
```

### Protein translation
```shell
ORF_finder \
-l 30 \
-n 10 \
-pc \
-log biopy/test_data/ORF_finder/ORF_finder.log \
-o biopy/test_data/ORF_finder/ \
biopy/test_data/ORF_finder/Ptrichocarpa_533_v4.1.cds.fa.gz
```

### RNA secondary structure prediction.
```python
from pybioinformatic import Nucleotide

# Generate random nucleic acid sequence.
random_nucl = Nucleotide.random_nucl(name='demo', length=[100, 150], bias=1.0)

# Secondary structure prediction
ss, mfe = random_nucl.predict_secondary_structure('biopy/test_data/RNA_structure/structure.ps')
print(ss, mfe, sep='\n')
```
```
>demo length=135
CAAAAAAAAACCATAAGCCGCCATGTCTCACATCGCAACCGGCTCAAGTAGAGTGCCCCTAATAATATGATCTTCGCTACAGAAGTTCCCCCCCCGCTGCCGGCTAGATGCGAACTCCACGCCTGGATGGCTCAG
...............((((((((.((......(((((.(.((((...((((.................((.((((......)))).))........)))).)))).).))))).......)).))).)))))...
-27.299999237060547
```
![image](test_data/RNA_structure/structure.png)

### Genotype consistency calculation.
```shell
gt_kit gs \
-i biopy/test_data/GT/GT.xls.gz \
-I biopy/test_data/GT/GT.xls.gz \
--database-compare \
-o biopy/test_data/GT/
```
![image](test_data/GT/Consistency.heatmap.png)

### Plot gene structure.
```shell
plot gene_structure \
-i biopy/test_data/gene_structure/Ptc.gff3.gz \
-o biopy/test_data/gene_structure/mRNA_structure.png
```
![image](test_data/gene_structure/mRNA_structure.png)

### Plot circos figure.
```shell
plot circos \
-c biopy/test_data/circos/Ptc_chr_len.txt \
-d biopy/test_data/circos/gene_density.txt \
-s biopy/test_data/circos/stat.txt \
-l biopy/test_data/circos/link.txt \
-o biopy/test_data/circos/circos.png
```
![image](test_data/circos/circos.png)

### Plot chromosome distribution.
```shell
plot chr_distribution \
-i biopy/test_data/chr_distribution/snp.ref.xls \
-l biopy/test_data/chr_distribution/chr_len.xls \
-w 100000 \
-n 4 \
-cmap RdYlGn \
-o biopy/test_data/chr_distribution/snp.distribution.png
```
![image](test_data/chr_distribution/snp.distribution.png)
