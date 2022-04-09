### Quick install and start ###
```
git clone https://github.com/zhangrengang/REPcluster
cd REPcluster

# install
conda -c bioconda kmer-db mcl
python3 setup.py install

# run an example
cd example_data
REPclust hifi.trf.fa -x 2	 # for tandem repeats
```

### Outputs ###
```
repclust.a2a.csv.jaccard	# Similariry matrix
repclust.network	# Network to import into Cytoscape
repclust.attr		# Attibutes to import into Cytoscape
repclust.mcl		# one Cluster per line
repclust.fa			# centered sequences for each cluster
```

