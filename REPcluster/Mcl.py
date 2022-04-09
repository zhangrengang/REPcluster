import sys

class MclGroup:
	def __init__(self, mcl):
		self.mcl = mcl
	def __iter__(self):
		return self._parse()
	def _parse(self):
		i = 0
		for line in open(self.mcl):
			i += 1
			id = 'C{}'.format(i)
			yield MclGroupRecord(line, id=id)

	def assign_cid(self, outLst, min_nodes=50):
		for grp in self:
			if len(grp) < min_nodes:
				continue
			for id in grp:
				lin = [id, grp.id]
				print('\t'.join(lin), file=outLst)

class MclGroupRecord:
	def __init__(self, line, id=None):
		self.id = id
		self.nodes = line.strip().split()
	def __len__(self):
		return len(self.nodes)
	def __iter__(self):
		return iter(self.nodes)

def main():
	MclGroup(sys.argv[1]).assign_cid(outLst=sys.stdout)

if __name__ == '__main__':
	main()
