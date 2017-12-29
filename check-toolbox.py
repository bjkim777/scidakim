import os
import sys
import subprocess
import argparse

## TOOL PATH ##
"""
def ToolPath(ver, tool, jar=None):
	ToolBox_path = 'Tools/Bio-tools/'

	if ver == 1.0:
		TOOL = { 'sickle': '1.33', 'bwa': '0.7.10',
		 'samtools':'1.1',  'picard': '1.68', 'bamqc': '2.1.3', 'gatk': '3.3',
		  'bamtools': '1.0.2', 'aneuploidy': '2.0.1', 'pigz': '2.3.3', 
		  'fastqValidator': '0.1.1a', 'bamValidator':'1.0.13', 'vcfValidator': '0.1.14'}

		return ToolBox_path + tool + '/' + TOOL[tool] + '/'
"""

def GetOpts():
	parser = argparse.ArgumentParser(description = 'Check ToolBox')
	parser.add_argument('-m', '--mount', type=str, required=True, metavar='<mount point>',  help= 'Toolbox mount point')
	parser.add_argument('-o', '--output', type=str, required=True,  metavar= '<output path>', help= 'Output path')
	argv = parser.parse_args()
	return(argv)

def Check(*files):
	for file in files:
		if os.path.isfile(file) and os.path.getsize(file) != 0:
			print (file + ' : OK')
			pass
		else:
			sys.stderr.write(
			"%s does not exist or gets size 0.\n" %file
			)	
			exit(1)

class MappingTool:
	def __init__(self, MOUNT):
		self.MOUNT = MOUNT

		## tools
		self.SICKLE = MOUNT + '/Tools/Bio-tools/sickle/1.33/sickle'
		self.BWA = MOUNT + '/Tools/Bio-tools/bwa/0.7.10/bwa'
		self.SAMTOOLS = MOUNT + '/Tools/Bio-tools/samtools/1.1/samtools'
		self.JAVA = MOUNT + '/Tools/Common/java/1.7.0_80/java'
		self.BAMVAL = MOUNT + '/Tools/Bio-tools/bamValidator/1.0.13/bam'
		self.FASTQVAL = MOUNT + '/Tools/Bio-tools/fastqValidator/0.1.1a/fastQValidator'
		self.VCFVAL = MOUNT + '/Tools/Bio-tools/vcfValidator/0.1.14/vcf-validator'
		self.BAMQC = MOUNT + '/Tools/Bio-tools/bamqc/2.1.3/qualimap'
		self.FASTQC = MOUNT + '/Tools/Bio-tools/fastqc/0.10.1/fastqc'
		self.JAR_GATK = MOUNT + '/Tools/Bio-tools/gatk/3.3/GenomeAnalysisTK.jar'
		self.JAR_MARKDUL = MOUNT + '/Tools/Bio-tools/picard/1.68/MarkDuplicates.jar'
		self.JAR_CALMET = MOUNT + '/Tools/Bio-tools/picard/1.68/CalculateHsMetrics.jar'

		## db
		self.REF = MOUNT + '/DB/Reference/fasta/hg19/1.0/human_ref.fa'
		self.DBSNP = MOUNT + '/DB/Reference/snv_indel/hg19/1.0/dbsnp_138.vcf'
		self.INDEL_VCF = MOUNT + '/DB/Reference/snv_indel/hg19/1.0/1000g_phase1_indels.vcf'
		self.MILLS = MOUNT + '/DB/Reference/snv_indel/hg19/1.0/mills_indel.vcf'
		self.P2 = MOUNT + '/DB/Region/panel/hg19/TS/1.0/P2.region'
	
	def GetFile(self):
		return {
		'name':'Check', 
		'fastq_r1':'Check_R1.fastq.gz', 
		'fastq_r2':'Check_R2.fastq.gz',
		}

	def WorkingTest(self, data_dir, output_path):
		print ('## Mapping Tool Test')

		## input data
		fastq_r1 = data_dir+'/'+self.GetFile()['fastq_r1']
		fastq_r2 = data_dir+'/'+self.GetFile()['fastq_r2']
		bait = data_dir + '/bait'

		## output data
		fastqc_result = output_path + '/' + self.GetFile()['name'] + '_R1_fastqc/fastqc_data.txt' 
		fastqval_result = output_path + '/' + self.GetFile()['name'] + '_fastq_val_result' 
		fastq_filterd_r1 = output_path + '/' + self.GetFile()['name'] + '_R1.fastq_filtered'
		fastq_filterd_r2 = output_path + '/' + self.GetFile()['name'] + '_R2.fastq_filtered'	
		sam = output_path + "/" + self.GetFile()['name'] + ".sam"
		bam = output_path + "/" + self.GetFile()['name'] + ".bam"
		sort_bam = output_path + "/" + self.GetFile()['name'] + ".sort.bam"
		rmdul_bam = output_path + "/" + self.GetFile()['name'] + ".sort_rmdul.bam"
		realign_target = output_path + "/" + self.GetFile()['name'] + ".intervals"
		realign_bam = output_path + "/" + self.GetFile()['name'] + ".sort_rmdul_realign.bam"
		table = output_path + "/" + self.GetFile()['name'] + ".table"
		recal_bam = output_path + "/" + self.GetFile()['name'] + ".sort_rmdul_realign_recal.bam"
		bamval_result = output_path + '/' + self.GetFile()['name'] + '_bamval_result'
		bamqc_result = output_path + '/' + 'genome_results.txt'
		metrix = output_path + '/metrix.txt'
		cov_met = output_path + '/cov_met.txt'

		## check tools and db
		print ('# Tools')
		Check(self.SICKLE, self.SAMTOOLS, self.BWA, self.JAVA, self.FASTQVAL, 
			self.BAMVAL, self.VCFVAL, self.FASTQC, self.BAMQC, self.JAR_GATK, self.JAR_MARKDUL, self.JAR_CALMET)
		print ('# DB' )
		Check(self.REF, self.DBSNP, self.INDEL_VCF, self.MILLS, self.P2)

		# FASTQC
		print ('# FASTQC')
		Check(fastq_r1, fastq_r2)
		subprocess.call([self.FASTQC, fastq_r1, fastq_r2, '-o', output_path])
		Check(fastqc_result)

		# FASTQVALIDATOR
		print ('# FASTQ Validator')
		p = subprocess.check_output([self.FASTQVAL, '--file', fastq_r1])
		with open(fastqval_result, 'w') as f:
			f.write(p)
		Check(fastqval_result)

		# trimming
		print ('# Trimming')
		subprocess.call([self.SICKLE, 'pe', '-f', fastq_r1, '-r', fastq_r2, '-o', fastq_filterd_r1, '-p', fastq_filterd_r2, 
			'-s', output_path + '/' + self.GetFile()['name'] + '_single.fastq_filtered', '-t', 'sanger', '-q', '20'])
		Check(fastq_filterd_r1, fastq_filterd_r2)

		# mapping
		print ('# Mapping')
		BAM_HEADER = "@RG\\tID:foo\\tLB:bar\\tPL:illumina\\tPU:illumina\\tSM:" + self.GetFile()['name']

		p = subprocess.Popen([self.BWA, 'mem', self.REF, fastq_filterd_r1, fastq_filterd_r2, '-R', BAM_HEADER], stdout = subprocess.PIPE)
		with open(sam, 'w') as f:
			for line in p.stdout.readlines():
				f.write(line)
		Check(sam)

		# sam to bam
		print ('# SAM to BAM')
		subprocess.call([self.SAMTOOLS, 'view', '-F', '0x800', '-bho', bam, sam])
		Check(bam)

		# sort
		print ('# Sorting')
		subprocess.call([self.SAMTOOLS, 'sort', '-T', output_path, '-o', sort_bam, bam])
		Check(sort_bam)

		# remove duplications
		print ('# Remove Duplications') 
		subprocess.call([self.JAVA, '-jar', self.JAR_MARKDUL, 'I=' + sort_bam, 'O=' + rmdul_bam, 'TMP_DIR=' + output_path, "M=" + output_path + "/" + self.GetFile()['name'] + ".met", "REMOVE_DUPLICATES=true"])
		Check(rmdul_bam)

		# index
		subprocess.call([self.SAMTOOLS, 'index', rmdul_bam])

		# create realign target
		print ('# Create target for realign') 
		subprocess.call([self.JAVA, '-jar', self.JAR_GATK, '-R', self.REF, '-T', 'RealignerTargetCreator', '-known', self.MILLS, '-known', self.INDEL_VCF, '-I', rmdul_bam, '-o', realign_target])
		Check(realign_target)

		# realign
		print ('# Realign') 
		subprocess.call([self.JAVA, '-jar', self.JAR_GATK, '-R', self.REF, '-T', 'IndelRealigner', '-I', rmdul_bam, '-targetIntervals', realign_target, '-known', self.MILLS, '-known', self.INDEL_VCF, '-o', realign_bam])
		Check(realign_bam)

		# recalibrate
		print ('# Recalibate') 
		subprocess.call([self.JAVA, '-jar', self.JAR_GATK, '-R', self.REF, '-T', 'BaseRecalibrator', '--knownSites', self.DBSNP, '--knownSites', self.MILLS, '--knownSites', self.INDEL_VCF, '-I', realign_bam, '-o', table])
		Check(table)

		# print read
		print ('# Print Reads') 
		subprocess.call([self.JAVA, '-jar', self.JAR_GATK, '-R', self.REF, '-T', 'PrintReads', '-BQSR', table, '-I', realign_bam, '-o', recal_bam])
		Check(recal_bam)

		subprocess.call([self.SAMTOOLS, 'index', recal_bam])

		# BAMQC
		print ('# BAMQC')
		Bin = os.getcwd()
		os.chdir(self.MOUNT + '/Tools/Bio-tools/bamqc/2.1.3/')
		subprocess.call([self.BAMQC, 'bamqc', '-bam', recal_bam,  '-outdir', output_path])
		Check(bamqc_result)
		os.chdir(Bin)

		# BAM Validator
		print ('# BAM Validator')
		p = subprocess.Popen([self.BAMVAL, 'validate', '--in', recal_bam], stderr=subprocess.PIPE)
		with open(bamval_result, 'w') as f:
			for line in p.stderr.readlines():
				f.write(line)
		Check(bamval_result)
		
		# metrix
		print ('# Matrix')
		subprocess.call([self.JAVA, '-jar', self.JAR_CALMET, 'R=' + self.REF,
			'BI=' + bait, 'TI=' + bait, 'N=' + self.GetFile()['name'], 'I=' + recal_bam, 'O=' + metrix, 
			'VALIDATION_STRINGENCY=LENIENT', 'PER_TARGET_COVERAGE=' + cov_met])
		Check(metrix, cov_met)

class SyntekaBioTool:
	def __init__(self, MOUNT):

		# tools
		self.ADISCAN = MOUNT + '/Tools/Gear/snv/1.0/adi_scan'
		self.CNVSCAN = MOUNT + '/Tools/Gear/cnv/1.0/cnv_scan'
		self.HAPLOSCAN = MOUNT + '/Tools/Gear/haplotype/1.0/haplo_scan'
		self.GHAPLOSCAN = MOUNT + '/Tools/Gear/haplotype/1.0/haplo_scan'
		self.ALLELECOUNT = MOUNT + '/Tools/Gear/haplotype/1.0/allele_count'
		self.RVR = MOUNT + '/Tools/Gear/indel/1.0/rvr'

		# db
		self.REF = MOUNT + '/DB/Reference/fasta/hg19/1.0/human_ref.fa'
		self.CNV = MOUNT + '/DB/Region/cnvscan/hg19/TS/1.0/cnv.region'

	def GetFile(self):
		return {
		'name':'Check',
		'bam': 'Check.sort_rmdul_realign_recal.bam',
		'hla': 'HLA',
		'rvr': 'DMET_d',
		'cnv': 'bam_list'
		}

	def WorkingTest(self, data_dir, output_path):
		print ('## SyntekaBio Tool Test')

		# input data
		bam = data_dir + "/" + self.GetFile()['bam']
		bam2 = output_path + "/" + self.GetFile()['bam']
		hla_path = data_dir + '/' + self.GetFile()['hla']
		rvr_d = data_dir + '/' + self.GetFile()['rvr']
		
		# output data
		fastq = output_path + "/" + self.GetFile()['name'] + '.fastq'
		filterd = output_path + '/' + self.GetFile()['bam'] + '.filtered'
		hla_gene = subprocess.Popen(['ls',hla_path], stdout=subprocess.PIPE).stdout.readlines()[0].strip()
		bam_list = output_path + '/' + hla_gene + '.list'
		hla_result = output_path + '/' + 'hla.report'
		rvr_result = output_path + '/' + 'rvr_result'
		adi_result = output_path + '/' + 'adi_result'
		cnv_list = output_path + '/' +  self.GetFile()['cnv']
		cnv_sex = output_path + '/' + 'sex'
		cnv_log = output_path + '/' + 'log'
		cnv_depth = output_path + '/' + 'depth'

		# check tools
		print ('# SyntekaBio Tools')
		Check(self.ADISCAN, self.CNVSCAN, self.HAPLOSCAN, self.GHAPLOSCAN, self.ALLELECOUNT, self.RVR)

		# ghaploscan
		print ('# Ghaplo_scan')
		p = subprocess.Popen([self.GHAPLOSCAN, '-a', bam, '-c', '1', '-s', '7000000', '-e', '8000000'], stdout=subprocess.PIPE)
		with open(fastq, 'w') as f:
			for line in p.stdout.readlines():
				f.write(line)
		Check(fastq)

		# allele count
		print ('# allele_count')
		subprocess.call([self.ALLELECOUNT, '-s', bam2, '-x', self.REF])
		Check(filterd)

		# haploscan
		print ('# haploscan')
		p = subprocess.Popen(['find', hla_path + '/' + hla_gene , '-name', '*.sorted.bam'], stdout=subprocess.PIPE)
		with open(bam_list, 'w') as f:
			for line in p.stdout.readlines():
				f.write(line)
		p = subprocess.Popen([self.HAPLOSCAN, bam_list, hla_gene, '0', '30', '1'], stdout=subprocess.PIPE)
		with open(hla_result, 'w') as f:
			for line in p.stdout.readlines():
				f.write(line)
		Check(bam_list, hla_result)

		# RVR
		print ('# rvr')
		p = subprocess.Popen([self.RVR, '-r', rvr_d, 'rrrn', '1'], stdout=subprocess.PIPE)
		with open(rvr_result, 'w') as f:
			for line in p.stdout.readlines():
				f.write(line)
		Check(rvr_result)

		# Adiscan
		print ('# adiscan')
		p = subprocess.Popen([self.ADISCAN, '-a', bam, '-r', self.REF, '-c', '2', '-s', '33141308', '-e', '33141553'], stdout=subprocess.PIPE)
		with open(adi_result, 'w') as f:
			for line in p.stdout.readlines():
				f.write(line)
		Check(adi_result)

		# CNVscan
		print ('# cvnscan')
		subprocess.call([self.CNVSCAN, '-v', cnv_list, '-l', self.CNV, '-o', output_path])
		Check(cnv_depth, cnv_log, cnv_sex)


### MAIN ###
def main (mount, output_path):
	if os.path.isdir(output_path):
		pass
	else:
		os.mkdir(output_path)

	check1 = MappingTool(mount)
	check2 = SyntekaBioTool(mount)
	check1.WorkingTest(mount+'/Check/input/', output_path)
	check2.WorkingTest(mount+'/Check/input/', output_path)

	#subprocess.call(['rm', '-rf', output_path])


### ACTION ### 
if __name__ == "__main__":
	argv = GetOpts()

	main(argv.mount, argv.output)
	sys.exit(0)