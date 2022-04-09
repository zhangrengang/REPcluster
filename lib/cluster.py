import sys,os
import copy
import argparse
import shutil
import json
import math
import multiprocessing
from collections import OrderedDict, Counter
from xopen import xopen as open
from Bio import SeqIO
from .small_tools import mkdirs, rmdirs, mk_ckp, check_ckp, test_s
from .RunCmdsMP import logger, run_cmd
from .multi_seqs import multi_seqs
from .matrix2list import matrix2list
from .Mcl import MclGroup
from .__version__ import version

NCPU = multiprocessing.cpu_count()
bindir = os.path.dirname(os.path.realpath(__file__))


def argparser():
	parser = argparse.ArgumentParser( 
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description='Cluster Repeat Sequences.',
		)
	# input
	group_in = parser.add_argument_group('Input', )
	group_in.add_argument('fasta', default=None, metavar='FILE', nargs='+', 
					help="Each sequence in a FASTA file is treated as a separate sample")
#	group_in.add_argument('-l', '-list', default=None, metavar='FILE', dest='lst',
#					help="Each line is a single FASTA file and treated as a separate sample")
	# output
	group_out = parser.add_argument_group('Output')
	group_out.add_argument('-pre', '-prefix', default='repclust', dest='prefix', metavar='STR',
					help="Prefix for output [default=%(default)s]")
	group_out.add_argument('-o', '-outdir', default='.', dest='outdir', metavar='DIR',
					help="Output directory [default=%(default)s]")
	group_out.add_argument('-tmpdir', default='tmp', type=str, metavar='DIR',
					help="Temporary directory [default=%(default)s]")
	# kmer
	group_kmer = parser.add_argument_group('Kmer matrix',)
	group_kmer.add_argument('-x', '-multiple', type=int, default=1, metavar='INT',
					dest='multiple',
					 help="Repeat sequences to cluster tandem repeat or circular sequences [default=%(default)s]")
	group_kmer.add_argument('-k', type=int, default=15, metavar='INT',
					 help="Length of kmer [default=%(default)s]")
	
	group_kmer.add_argument('-m', "-measure", default='jaccard', dest='measure',
					choices=['jaccard', 'min', 'max', 'cosine',],
					help="The similarity measure to be calculated. [default=%(default)s]")

 
	# cluster
	group_clst = parser.add_argument_group('Cluster', )
	group_clst.add_argument('-c', '-min_similarity', type=float, default=0.2, metavar='FLOAT', 
					dest='min_similarity',
					help="Minimum similarity to cluster [default=%(default)s]")
					
	group_clst.add_argument('-I', '-inflation', type=float, default=2.0, metavar='FLOAT',
					dest='inflation',
					help="Inflation for MCL (varying this parameter affects granularity) [default=%(default)s]")

	# others
	group_other = parser.add_argument_group('Other options')
	group_other.add_argument('-p', '-ncpu', type=int, default=NCPU, metavar='INT', dest='ncpu',
					 help="Maximum number of processors to use [default=%(default)s]")
	group_other.add_argument('-cleanup', action="store_true", default=False,
					help="Remove the temporary directory [default=%(default)s]")	
	group_other.add_argument('-overwrite', action="store_true", default=False,
					help="Overwrite even if check point files existed [default=%(default)s]")
	group_other.add_argument('-v', '-version', action='version', version=version)
	return parser

def makeArgparse():
	parser = argparser()
	args = parser.parse_args()
	if args.prefix is not None:
		args.prefix = args.prefix.replace('/', '_')
	return args

class Pipeline:
	def __init__(self, **kargs):
		self.__dict__.update(**kargs)
		
	def run(self):
		# mkdir
		self.outdir = os.path.realpath(self.outdir)
		self.tmpdir = os.path.realpath(self.tmpdir)
		mkdirs(self.outdir)
		mkdirs(self.tmpdir)
		self.outdir += '/'
		self.tmpdir += '/'
		self._outdir = self.outdir
		if self.prefix is not None:
			self.outdir = self.outdir + self.prefix
			self.tmpdir = self.tmpdir + self.prefix

		logger.info('Build distance matrix..')
		
		opts = '-k {}'.format(self.k)
		fasta = self.tmpdir + '.multi.fa'
		with open(fasta, 'w') as fout:
			d_seqs = multi_seqs(seqfiles=self.fasta, outfile=fout, 
						fold=self.multiple, min_length=self.k*self.multiple)
		opts += ' -multisample-fasta'
		input = self.tmpdir + '.list'
		cmd = 'realpath {} > {}'.format(fasta, input)
		run_cmd(cmd, log=True)
			
		# use kmer-db for marix
		db = self.tmpdir + '.kmer.db'
		matrix = self.outdir + '.a2a.csv'
		cmd = '''kmer-db build {opts} {input} {db} && \
kmer-db all2all {db} {matrix} && \
kmer-db distance {measure} -phylip-out {matrix}'''.format(
			opts=opts, measure=self.measure, input=input, db=db, matrix=matrix)
		run_cmd(cmd, log=True)
		
		dist = matrix + '.' + self.measure
		network = self.outdir + '.network'
		with open(network, 'w') as fout:
			matrix2list(dist, fout, cutoff=self.min_similarity, phylip=True)
		
		# cluster by mcl
		logger.info('Cluster..')
		cluster = self.outdir + '.mcl'
		cmd = 'mcl {input} --abc -I {inflation} -o {output} -te {ncpu}'.format(
			inflation=self.inflation, input=network, output=cluster, ncpu=self.ncpu)
		run_cmd(cmd, log=True)
		
		# attribution of CID
		attr = self.outdir + '.attr'
		mcl = MclGroup(cluster, prefix='CL')
		with open(attr, 'w') as fout:
			mcl.assign_cid(fout, min_nodes=0, )
		logger.info('Import `{}` and `{}` into Cytoscape for visualization'.format(network, attr))

		# output seq
		outseq = self.outdir + '.fa'
		logger.info('Output centred sequences: `{}`'.format(outseq))
		
		with open(outseq, 'w') as fout:
			mcl.generate_seqs(fout, network, d_seqs)
		
def main(args=makeArgparse(), opts=None):
	if args is None and isinstance(opts, str):
		parser = argparser()
		args = parser.parse_args(opts.split())
	logger.info('Command: {}'.format(' '.join(sys.argv)))
	logger.info('Version: {}'.format(version))
	logger.info('Arguments: {}'.format(args.__dict__))
	pipeline = Pipeline(**args.__dict__)
	return pipeline.run()


if __name__ == '__main__':
	main()