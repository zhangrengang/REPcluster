import sys
import itertools

class MclGroup:
	def __init__(self, mcl, prefix='C'):
		self.mcl = mcl
		self.prefix = prefix
	def __iter__(self):
		return self._parse()
	def _parse(self):
		i = 0
		for line in open(self.mcl):
			i += 1
			id = '{}{}'.format(self.prefix, i)
			yield MclGroupRecord(line, id=id)

	def assign_cid(self, outLst, min_nodes=50):
		for grp in self:
			if len(grp) < min_nodes:
				continue
			for id in grp:
				lin = [id, grp.id]
				print('\t'.join(lin), file=outLst)

	def assign_center(self, abc):
		d_weight = MclInput(abc).to_dict()
		for grp in self:
			d_sum = {}
			for node1, node2 in grp.comb:
				try: weight = d_weight[(node1, node2)]
				except KeyError: continue
				for node in (node1, node2):
					try: d_sum[node] += weight
					except KeyError: d_sum[node] = weight
			if len(d_sum) > 1:
				grp.center = max(d_sum.items(), key=lambda x:x[1])[0]
			else:
				grp.center = grp[0]
			yield grp
		
	def generate_seqs(self, fout, abc, d_seqs):
		for grp in self.assign_center(abc):
			seq = d_seqs[grp.center]
			desc = 'N={};L={};center={}'.format(len(grp), len(seq), grp.center)
			print('>{} {}\n{}'.format(grp.id, desc, seq), file=fout)
			
class MclGroupRecord:
	def __init__(self, line, id=None):
		self.id = id
		self.nodes = line.strip().split()
	def __len__(self):
		return len(self.nodes)
	def __getitem__(self, idx):
		return self.nodes[idx]
	def __iter__(self):
		return iter(self.nodes)
	@property
	def comb(self, ):
		return itertools.combinations(self, 2)

class MclInput:
	def __init__(self, abc):
		self.abc = abc
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in open(self.abc):
			yield MclInputLine(line)
	def to_dict(self):
		d = PairDict()
		for line in self:
			d[line.key] = line.weight
		return d
class MclInputLine:
	def __init__(self, line):
		node1, node2, weight = line.strip().split() 
		self.node1, self.node2 = node1, node2
		self.weight = float(weight)
	@property
	def key(self):
		return (self.node1, self.node2)
	def in_dict(self, d):
		return (self.node1, self.node2) in d or (self.node2, self.node1) in d

class PairDict(dict):
	def __init__(self, *args, **kargs):
		super().__init__(*args, **kargs)
	def __setitem__(self, key, value, dict_setitem=dict.__setitem__):
		key = self.format_key(key)
		return dict_setitem(self, key, value)
	def __getitem__(self, key, dict_getitem=dict.__getitem__):
		key = self.format_key(key)
		return dict_getitem(self, key)
	def __delitem__(self, key, dict_delitem=dict.__delitem__):
		key = self.format_key(key)
		dict_delitem(self, key)
	
	def format_key(self, key):
		return tuple(sorted(key))
def main():
	MclGroup(sys.argv[1]).assign_cid(outLst=sys.stdout)

if __name__ == '__main__':
	main()
