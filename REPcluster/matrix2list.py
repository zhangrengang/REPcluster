import sys

def matrix2list(inMat, outLst, cutoff=0.2, phylip=False):
	i = 0
	last_ids = []
	for line in open(inMat):
		i += 1
		if i == 1:
			continue
		if phylip:
			temp = line.rstrip().split()
		else:
			temp = line.rstrip('\n,').split(',')
		id = temp[0].split()[0]
		s = len(temp) - (i-2) 
		values = map(float, temp[s:])
		for last_id, value in zip(last_ids, values):
			if value < cutoff:
				continue
			line = [last_id, id, value]
			line = map(str, line)
			print('\t'.join(line), file=outLst)
		last_ids += [id]

def main():
	try: cutoff = float(sys.argv[2])
	except IndexError: cutoff = 0.2
	matrix2list(inMat=sys.argv[1], outLst=sys.stdout, cutoff=cutoff)

if __name__ == '__main__':
	main()
