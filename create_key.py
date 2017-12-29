import sys


def MkKey(label , chro , start , end):
	with open(label + ".key" , 'w') as key :
		for num in range(int(start) , int(end)+1):
			key.write(chro + "_" + str(num) + "\n")
	print (label)
	return key



### MAIN ###
def main (boundary):
	with open(boundary,'r') as f:
		for line in f.readlines():	
			col=line.strip().split('\t')
			label=col[0]
			chro=col[1]
			start=col[2]
			end=col[3]

			MkKey(label, chro, start, end)
			

### ACTION ### 
if __name__ == "__main__":
	if len(sys.argv) == 2:
		 boundary= sys.argv[1]
	else:
		sys.stderr.write(
"Usage : "+sys.argv[0]+" [boundary]\n"
)
		sys.exit(1)

	main( boundary)
	sys.exit(0)