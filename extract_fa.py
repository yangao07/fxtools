from pyfastx import Fasta
import sys

def get_names(name_fn):
	d = dict()
	with open(name_fn) as in_fp:
		for line in in_fp:
			d[line.rstrip()] = 1
	# print(d)
	return d

def extract_fa(names, in_fa):
	for name in in_fa.keys():
		# print(name)
		if name in names:
			print('>{}\n{}'.format(name, in_fa[name][:].seq.upper()))

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print('{} in.name in.fa > out.fa'.format(sys.argv[0]))
		sys.exit(1)

	names = get_names(sys.argv[1])
	fa = Fasta(sys.argv[2])

	extract_fa(names, fa)