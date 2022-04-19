REPcluster aim to cluster repeat sequences that have similar contents but distinct structures, such as (TTTAGGG)m vs (TTTAGGG)n (tandem repeats), A-B-C vs A-C (some TE sequences).

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

### Usage ###
```
usage: REPclust [-h] [-pre STR] [-o DIR] [-tmpdir DIR] [-x INT] [-k INT]
               [-m {jaccard,min,max,cosine}] [-c FLOAT] [-I FLOAT] [-p INT]
               [-cleanup] [-overwrite] [-v]
               FILE [FILE ...]

Cluster Repeat Sequences.

optional arguments:
  -h, --help            show this help message and exit

Input:
  FILE                  Each sequence in a FASTA file is treated as a separate
                        sample

Output:
  -pre STR, -prefix STR
                        Prefix for output [default=repclust]
  -o DIR, -outdir DIR   Output directory [default=.]
  -tmpdir DIR           Temporary directory [default=tmp]

Kmer matrix:
  -x INT, -multiple INT
                        Repeat sequences to cluster tandem repeat or circular
                        sequences [default=1]
  -k INT                Length of kmer [default=15]
  -m {jaccard,min,max,cosine}, -measure {jaccard,min,max,cosine}
                        The similarity measure to be calculated.
                        [default=jaccard]

Cluster:
  -c FLOAT, -min_similarity FLOAT
                        Minimum similarity to cluster [default=0.2]
  -I FLOAT, -inflation FLOAT
                        Inflation for MCL (varying this parameter affects
                        granularity) [default=2.0]

Other options:
  -p INT, -ncpu INT     Maximum number of processors to use [default=32]
  -cleanup              Remove the temporary directory [default=False]
  -overwrite            Overwrite even if check point files existed
                        [default=False]
  -v, -version          show program's version number and exit
```
