### Quick install and start ###
```
git clone https://github.com/zhangrengang/REPcluster
cd REPcluster

# install
conda -c bioconda kmer-db mcl
python3 setup.py install

# run an example
cd example_data
REPclust -f hifi.high.trf -m min -k 13 -c 0.1 -x 2
```
