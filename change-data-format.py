import sys
import argparse
import re

def GetOpts():
	parser = argparse.ArgumentParser(description = '')
	parser.add_argument('-d', '--data', type=str, required=True,  metavar= '<data path>', help= 'data path')
	argv = parser.parse_args()
	return(argv)

def changeData(data):
	with open(data, 'r') as f:
		PATTEN = ""
		for line in f:
			if re.match("-",line):
				size = line.strip().split()[4]
				date = "/".join(line.strip().split()[5:8])
				data = line.strip().split()[-1]

				if PATTEN:
					print "{path}/{file}\t{size}\t{date}".format(path=PATTEN, file=data, size=size, date=date)
				else:
					print "{file}\t{size}\t{date}".format(file=data, size=size, date=date)

			elif len(line.strip().split())==1 and re.search(":$", line):
				PATTEN = line.strip().split(':')[0]

### MAIN ###
def main (data):
	changeData(data)
			

### ACTION ### 
if __name__ == "__main__":
	argv = GetOpts()

	main( argv.data )
	sys.exit(0)