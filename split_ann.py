import sys
import subprocess
import tempfile

###### TOOL ######
BEDTOOLS="/data/tools/bedtools2/bin/bedtools" # v2.26

class Region:
	def __init__ (self, bed_chr,bed_start,bed_end,bed_name,file,ref):
		self.chr=bed_chr
		self.start=bed_start
		self.end=bed_end
		self.label=bed_name
		self.file=file
		self.ref=ref

	### CHOICE ###
	def Action(self):
		# extract first position
		head = subprocess.Popen(['head','-n',self.start,self.file] , stdout=subprocess.PIPE)
		tail = subprocess.Popen(['tail','-1'] , stdin=head.stdout,stdout=subprocess.PIPE)
		head.stdout.close()
		first_line = tail.communicate()[0].strip()
		chro=first_line.split()[2]
		first_pos=first_line.split()[3]

		# decide to exist start position
		awk = subprocess.Popen(['awk','{print $4}',self.file],stdout=subprocess.PIPE)
		isStart = subprocess.call(['grep','-w',self.start],stdin=awk.stdout,stdout=subprocess.PIPE)
		awk.stdout.close()

		# decide to exist end position
		isEnd=subprocess.call(['grep','-w',self.end,self.file] ,stdout=subprocess.PIPE)

		if self.chr == chro :
			if self.start == first_pos :
				if isEnd == 0 :											# full postions
					return self.SplitAnn(self.start,self.end)
				else :													# half(tail)
					# extract position in end line
					end_line = subprocess.Popen(['tail','-1',self.file] ,stdout=subprocess.PIPE)
					end_pos=end_line.stdout.readlines()[0].strip().split()[3]

					return self.AddRef(end=end_pos)
			else :
				if isEnd != 0 :												# no positions
					return self.AddRef()
				else :
					if isStart != 0 :										# half(head)
						return self.AddRef(start=first_pos)
					else : 												# split into middle
						# extract position in first line
						sub_line = subprocess.Popen(['head','-1',self.file] ,stdout=subprocess.PIPE)
						file_first_pos = sub_line.stdout.readlines()[0].strip().split()[3]

						# find line number
						middle_start_pos = int(self.start)-int(file_first_pos)+1
						middle_end_pos = int(self.end)-int(file_first_pos)+1

						return self.SplitAnn(str(middle_start_pos) , str(middle_end_pos))
		else:
			#sys.stderr.write('Chromosome is not matching.\n')
			pass

	### CASE FANCTION ###
	def SplitAnn(self,start,end):
		region=start+","+end+"p"
		p1=subprocess.Popen(['sed','-n',region,self.file] , stdout=subprocess.PIPE)

		with open(self.label+".ann",'w') as ann:
			for line in p1.stdout.readlines():
				ann.write(line)
		return ann


	def AddRef(self,start=None,end=None):
		if start :
			num=int(start)-1
			ann=self.MkAnn(self.chr , self.start , num)
			with open(ann.name, 'a') as add_ann:
				region = str(int(self.end ) - int(start) +1)
				p1 = subprocess.Popen(['head','-n',region,self.file] , stdout=subprocess.PIPE)
				for line in p1.stdout.readlines() :
					add_ann.write(line)
			return add_ann

		elif end :		# not complete
			sys.stderr.write("Function of no end position in ann file  is not complete.")
			sys.stderr.write("boundary : "+self.chr+" "+self,start+" "+self.end)
			sys.exit(1)
			"""
			num=int(end)+1
			ann=self.MkAnn(self.chr,str(num),self.end)
			with opne(ann, 'a') as add_ann :
				region = str(int(end) - int(self.start) +1)
				p1 = subprocess.Popen(['tail','-n',region,self.file] , stdout=subprocess.PIPE)
				for line in p1.stdout.readlines() :
					add_ann.write(line)
			"""
			return add_ann
		else:
			ann=MkAnn(self.chr,self.start,self.end)
			return ann


	### ACTIVITY FUNCTION ###
	def MkAnn(self,chro,start,end):
		bed = self.MkBed(chro,start,end)
		getfasta = subprocess.Popen([ BEDTOOLS , 'getfasta' , '-fi' , self.ref ,'-bed' , bed.name, '-tab' ] , stdout=subprocess.PIPE)
		seq = getfasta.stdout.readlines()[0].strip().split()[1]
		
		with open(self.label+".ann", 'w') as ann:
			plus=int(start)
			for seq in list(seq) : 
				if seq == 'N' :
					ann.write(self.chr + "_" + str(plus) + " Intergenic "+chro+" "+str(plus)+" "+"none_"+str(plus)+" "+seq+" A\n")
				elif seq == "A" :
					ann.write(self.chr + "_" + str(plus) + " Intergenic "+chro+" "+str(plus)+" "+"none_"+str(plus)+" "+seq+" T\n")
				elif seq == 'T' :
					ann.write(self.chr + "_" + str(plus) + " Intergenic "+chro+" "+str(plus)+" "+"none_"+str(plus)+" "+seq+" A\n")
				elif seq == 'G' :
					ann.write(self.chr + "_" + str(plus) + " Intergenic "+chro+" "+str(plus)+" "+"none_"+str(plus)+" "+seq+" C\n")
				elif seq == 'C' :
					ann.write(self.chr + "_" + str(plus) + " Intergenic "+chro+" "+str(plus)+" "+"none_"+str(plus)+" "+seq+" G\n")
				plus+=1
		
		# remove temporary file(bed)
		subprocess.call(['rm','-f',bed.name])

		return ann
		

	def MkBed(self,chro,start,end):
		temp = tempfile.NamedTemporaryFile(delete=False,bufsize=1)
		minus=int(start)-1
		temp.write(str(chro)+"\t"+str(minus)+"\t"+str(end)+"\n")
		return temp


### MAIN ###
def main (boundary,ann,ref):
	with open(boundary,'r') as f:
		for line in f.readlines():	
			col=line.strip().split('\t')
			label=col[0]
			chro=col[1]
			start=col[2]
			end=col[3]

			data=Region(chro,start,end,label,ann,ref)
			data.Action()
			

### ACTION ### 
if __name__ == "__main__":
	if len(sys.argv) == 4:
		 boundary= sys.argv[1]
		 ann= sys.argv[2]
		 ref= sys.argv[3]
	else:
		sys.stderr.write(
"Usage : "+sys.argv[0]+" [boundary] [ann] [ref]\n"
"""
If Ann file has gap(chromosome position),  do not use this script.
This script needs bedtools v2.26.

Boundary format
-------------------------
[label]	[chr]	[start pos]	[end pos]
ex) 1_1_1	1	1	4599925

Ann format
-------------------------
[chr_pos] [function] [chr] [position] [rs id] [ref seq] [alt seq]
ex) 1_5 Intergenic 1 5 none_5 N A

Ref is Fasta.

"""
)
		sys.exit(1)

	main( boundary, ann,ref)
	sys.exit(0)