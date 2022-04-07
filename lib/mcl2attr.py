import sys

def assign_cid(inMcl, outLst, min_nodes=50):
	i = 0
	for line in open(inMcl):
		i += 1
		cid = 'C{}'.format(i)
		ids = line.strip().split()
		if len(ids) < min_nodes:
			continue
#		print >> sys.stderr, cid, len(ids)
		for id in ids:
			lin = [id, cid]
			print('\t'.join(lin), file=outLst)

def main():
	assign_cid(inMcl=sys.argv[1], outLst=sys.stdout)

if __name__ == '__main__':
	main()
