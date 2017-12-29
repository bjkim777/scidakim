import os
import sys
import subprocess
import tempfile
import argparse
import sqlite3

EXOME_ANN_DB = 'exome-annotation.db'
AMINOACID = {
	"C":"Cys", "S":"Ser", "T":"Thr", "P":"Pro", "A":"Ala", 
	"G":"Gly", "N":"Asn", "D":"Asp", "E":"Glu", "Q":"Gln",
	"H":"His", "R":"Arg", "K":"Lys", "M":"Met", "I":"Ile",
	"L":"Leu", "V":"Val", "F":"Phe", "Y":"Tyr", "W":"Trp", "X":"Stop"}

def GetOpts():
	parser = argparse.ArgumentParser(description = '')
	parser.add_argument('-a', '--arg1', type=str, default='.', metavar='<arg1>',  help= 'description')
	parser.add_argument('-b', '--arg2', type=str, required=True,  metavar= '<arg2>', help= 'description')
	argv = parser.parse_args()
	return(argv)

"""
class ANN:
	def __init__(self, output, adi):
		self.output = output
		self.adi = adi
"""

def AnnAnnovar(ADI):
	ANNOVAR = '/tools/Bio-tools/annovar/table_annovar.pl'
	ANNOVAR_DB = '/tools/Bio-tools/annovar/humandb'

	#ANNOVAR_INPUT_FILE = tempfile.NamedTemporaryFile(delete=False, bufsize=1)
	ANNOVAR_INPUT_FILE =  "1." + str(ADI).split('/')[-1].split('.')[0]  + ".input"
	ANN_ANNOVAR = str(ANNOVAR_INPUT_FILE) + ".hg19_multianno.txt"
	ANN_ANNOVAR_AFTER = "1.2." + str(ADI).split('.')[0] + ".result"

	DAMAGE = {"D":"probably_damaging", "P":"possibly_damaging", "B":"benign", "NA":"-", ".":"-"}
	GYS2IHT = {"ABCA12":"recessive", "ACTG1":"dominant", "ADA":"recessive", "AFG3L2":"recessive", "AIRE":"dominant", "ALDH18A1":"recessive", "ALOX12B":"recessive", "ALOXE3":"recessive", "ANKH":"dominant", "ANO10":"recessive", "APOE":"Mixed", "AQP2":"Mixed", "ARFGEF2":"recessive", "ARHGEF10":"dominant", "ATL1":"dominant", "ATOH7":"recessive", "ATP6V0A4":"recessive", "BCO1":"dominant", "BEST1":"recessive", "BIN1":"recessive",
		"CA2":"recessive", "CASR":"dominant", "CAV3":"recessive", "CDH15":"dominant", "CDH23":"recessive", "CDK6":"recessive", "CENPE":"recessive", "CENPJ":"recessive", "CEP152":"recessive", "CHCHD10":"dominant", "CIB2":"recessive", "CLCN1":"Mixed", "CLCN7":"Mixed", "CLDN14":"recessive", "CLN6":"recessive", "COCH":"dominant", "COL11A2":"Mixed", "COL2A1":"dominant", "COL4A3":"dominant", "COL4A4":"recessive",
		"COL6A1":"dominant", "COL6A2":"Mixed", "COL7A1":"Mixed", "CPT1C":"dominant", "CRADD":"recessive", "CRYBB3":"Mixed", "CRYM":"dominant", "CYBA":"recessive", "CYCS":"dominant", "CYP2U1":"recessive", "CYP4F22":"recessive", "CYP7B1":"recessive", "DCAF8":"dominant", "DDHD2":"recessive", "DEAF1":"dominant", "DFNB59":"recessive", "DIABLO":"dominant", "DLL3":"recessive", "DNAJC5":"dominant", "DNMT1":"dominant",
		"DOCK8":"recessive", "DSPP":"dominant", "DYNC1H1":"dominant", "EEF1A2":"dominant", "EFEMP2":"recessive", "EGR2":"Mixed", "EIF4G1":"dominant", "ELANE":"dominant", "ELMOD3":"recessive", "ENPP1":"recessive", "ESPN":"dominant", "FA2H":"recessive", "FAS":"recessive", "FBLN5":"recessive", "FBN1":"dominant", "FBXO7":"recessive", "FGF23":"dominant", "FUS":"recessive", "FYCO1":"recessive", "G6PC3":"recessive",
		"GBA2":"recessive", "GCH1":"recessive", "GDAP1":"Mixed", "GDF3":"dominant", "GDF6":"dominant", "GFI1":"dominant", "GIGYF2":"dominant", "GIPC3":"recessive", "GJA1":"recessive", "GJB2":"Mixed", "GJB3":"Mixed", "GJB6":"dominant", "GLRA1":"Mixed", "GNA11":"dominant", "GNAT1":"dominant", "GNE":"recessive", "GPT2":"recessive", "GRIN1":"dominant", "GRIN2B":"dominant", "GRXCR1":"recessive",
		"HCN4":"dominant", "HERC2":"recessive", "HINT1":"recessive", "HPCA":"recessive", "HPGD":"recessive", "HTRA1":"recessive", "HTRA2":"dominant", "IGHMBP2":"recessive", "IGLL1":"recessive", "IL7R":"recessive", "JAGN1":"recessive", "JAK3":"recessive", "KARS":"recessive", "KCNQ4":"dominant", "KIAA0196":"dominant", "KIAA1033":"recessive", "KIF1A":"Mixed", "KIF1C":"recessive", "KIRREL3":"dominant", "KLHL3":"Mixed",
		"KRT14":"recessive", "KRT5":"recessive", "KRT74":"dominant", "LFNG":"recessive", "LGI1":"dominant", "LHFPL5":"recessive", "LIM2":"recessive", "LIPH":"recessive", "LMNA":"Mixed", "LPAR6":"recessive", "LRP5":"Mixed", "LRP6":"dominant", "LRRK2":"dominant", "LRTOMT":"recessive", "LTBP4":"recessive", "MAB21L2":"recessive", "MAN1B1":"recessive", "MAT1A":"Mixed", "MEFV":"dominant", "MESP2":"recessive",
		"MMP13":"dominant", "MPZ":"dominant", "MSRB3":"recessive", "MTMR14":"dominant", "MTPAP":"recessive", "MYH14":"dominant", "MYH7":"Mixed", "MYH9":"dominant", "MYO15A":"recessive", "MYO1A":"dominant", "MYO6":"recessive", "MYO7A":"Mixed", "NCF1":"recessive", "NCF2":"recessive", "NCF4":"recessive", "NDST1":"recessive", "NHP2":"recessive", "NIPA1":"dominant", "NIPAL4":"recessive", "NOP10":"recessive",
		"NOTCH3":"dominant", "NR3C2":"dominant", "NRL":"recessive", "NSUN2":"recessive", "OPA3":"dominant", "OTOF":"recessive", "P2RX2":"dominant", "PACS1":"dominant", "PARK2":"recessive", "PARK7":"recessive", "PARN":"recessive", "PCDH15":"recessive", "PDE6B":"dominant", "PEO1":"dominant", "PHC1":"recessive", "PINK1":"recessive", "PKHD1":"recessive", "PLEKHG5":"recessive", "PMP22":"Mixed", "PNPLA6":"recessive",
		"PNPT1":"recessive", "POLG2":"dominant", "POLG":"Mixed", "POU4F3":"dominant", "PPP2R1A":"dominant", "PRIMPOL":"dominant", "PROC":"Mixed", "PROS1":"Mixed", "PURA":"dominant", "PYCR1":"recessive", "RDH5":"recessive", "RDX":"recessive", "REEP1":"dominant", "REEP2":"Mixed", "RHO":"Mixed", "RNF170":"dominant", "ROR2":"recessive", "RTEL1":"Mixed", "RTN2":"dominant", "RYR1":"recessive",
		"SCN5A":"recessive", "SCNN1A":"recessive", "SCNN1B":"recessive", "SETX":"recessive", "SLC16A1":"dominant", "SLC17A8":"dominant", "SLC25A38":"recessive", "SLC25A4":"dominant", "SLC26A4":"recessive", "SLC2A1":"recessive", "SLC33A1":"dominant", "SLC39A5":"dominant", "SLC4A1":"dominant", "SLC6A17":"recessive", "SLC6A5":"Mixed", "SLCO2A1":"recessive", "SMARCA4":"dominant", "SMARCB1":"dominant", "SNCA":"dominant", "SNX10":"recessive",
		"SOD1":"recessive", "SOST":"dominant", "SOX11":"dominant", "SPG7":"recessive", "SPTA1":"recessive", "SPTB":"dominant", "ST14":"recessive", "ST3GAL3":"recessive", "STAT3":"dominant", "STIL":"recessive", "STUB1":"recessive", "SYNGAP1":"dominant", "SYT14":"recessive", "TAF2":"recessive", "TBC1D24":"Mixed", "TCIRG1":"recessive", "TDP1":"recessive", "TECR":"recessive", "TECTA":"dominant", "TERT":"Mixed",
		"TFG":"recessive", "TGM1":"recessive", "THRB":"dominant", "TH":"recessive", "TINF2":"dominant", "TMC1":"Mixed", "TMEM43":"dominant", "TMIE":"recessive", "TMPRSS3":"recessive", "TNC":"dominant", "TNFRSF11A":"recessive", "TNFRSF1A":"dominant", "TNFSF11":"recessive", "TOR1A":"dominant", "TPP1":"recessive", "TREX1":"dominant", "TTI2":"recessive", "TUBB1":"dominant", "TUBB4A":"dominant", "UCHL1":"dominant",
		"VPS37A":"recessive", "WDR62":"recessive", "WFS1":"dominant", "WNT5A":"dominant", "WRAP53":"recessive", "WWOX":"recessive", "ZFYVE27":"dominant", "ZNF335":"recessive", "ZNF592":"recessive", "ZNF644":"dominant"}


	ANNOVAR_INPUT = subprocess.Popen(
		['awk', '{split($2,pos,"_"); \
		if (length($5)!=1) \
		{print pos[1], pos[2], pos[2]+length($5)-1, $5, $6} \
		else  \
		{print pos[1], pos[2], pos[2], $5, $6} }', ADI], stdout=subprocess.PIPE)

	with open(ANNOVAR_INPUT_FILE,'w') as f:
		for line in ANNOVAR_INPUT.stdout.readlines():
			f.write(line)

	subprocess.call(['perl', str(ANNOVAR), 
		str(ANNOVAR_INPUT_FILE), 
		str(ANNOVAR_DB), 
		'--buildver', 'hg19', 
		'-remove', 
		'-protocol', 'refGene,avsnp147,ljb26_all', 
		'-operation', 'g,f,f', 
		'-nastring', 'NA', 
		'-otherinfo'])

	REMOVE_HEAD=subprocess.Popen(['sed','1d',ANN_ANNOVAR], stdout=subprocess.PIPE)

	for line in REMOVE_HEAD.stdout.readlines():
		COL = line.strip().split('\t')

		# COLUMN
		CHR = COL[0]
		START = COL[1]
		REF = COL[3]
		ALT = COL[4]
		FUNC = COL[5]
		GENE = COL[6]
		GENE_DETAIL = COL[7]
		EXONICFUNC = COL[8].replace(' ', '_')
		AACHANGE = COL[9]
		DBSNP147 = COL[10]
		SIFT_SCORE = COL[11]
		HVAR_SCORE = COL[15]
		HVAR_PRED = COL[16]
		TRN = 0

		# Change easy word
		if HVAR_SCORE == '.':
			HVAR_SCORE = 'NA'

		# Variant Class
		if len(REF) == len(ALT):
			VAR_CLA = 'SNV'
		else:
			VAR_CLA = 'INDEL'

		# Dominant Recessive
		if len(GENE.split(',')) == 1:
			if GENE in GYS2IHT:
				IHT = GYS2IHT[str(GENE)]
			else :
				IHT = 'other'
		else:
			IHT = 'other'

		# Splicing OX
		if FUNC.find('splicing')!=-1:
			SPL = 1
		else:
			SPL = 0

		# Extract Transcript ID
		if AACHANGE != 'NA' and AACHANGE != 'UNKNOWN':
			TRANS = AACHANGE.split(',')[0].split(':')[1]

		# Print line 
		if EXONICFUNC != 'unknown' and EXONICFUNC != 'NA':
			if len(AACHANGE.split(',')[0].split(':')) == 5: # amino acid change OX
				AMINO_INFO = AACHANGE.split(',')[0].split(':')[-1].split('.')[-1]

				if AMINO_INFO[0] == 'X' or AMINO_INFO[-1] == 'X':
					TRN = 1

				RESULT="{chr}_{pos}\t{ref}\t{alt}\t{var_cla}\t{dbsnp147}\t{iht}\t{func}\t{gene}\t{exonicfunc}\t{trans_id}\tp.{amino1}{region}{amino2}\t{aachange}\t{spl}\t{trn}\t{hvar_score}\t{sift}\t{damage}\n".format(
					chr=CHR, pos=START, ref=REF, alt=ALT,
					var_cla=VAR_CLA, dbsnp147=DBSNP147, iht=IHT, func=FUNC,
					gene=GENE, exonicfunc=EXONICFUNC, trans_id=TRANS, amino1=AMINOACID[AMINO_INFO[0]],
					region=AMINO_INFO[1:-1], amino2=AMINOACID[AMINO_INFO[-1]], aachange=AACHANGE, spl=SPL,
					trn=TRN, hvar_score=HVAR_SCORE, sift=SIFT_SCORE, damage=DAMAGE[str(HVAR_PRED)])

				with open(ANN_ANNOVAR_AFTER, 'a') as f:
					f.write(RESULT)
			else:
				RESULT="{chr}_{pos}\t{ref}\t{alt}\t{var_cla}\t{dbsnp147}\t{iht}\t{func}\t{gene}\t{exonicfunc}\t{trans_id}\t-\t{aachange}\t{spl}\t{trn}\t{hvar_score}\t{sift}\t{damage}\n".format(
					chr=CHR, pos=START, ref=REF, alt=ALT,
					var_cla=VAR_CLA, dbsnp147=DBSNP147, iht=IHT, func=FUNC,
					gene=GENE, exonicfunc=EXONICFUNC, trans_id=TRANS, aachange=AACHANGE,
					spl=SPL, trn=TRN, hvar_score=HVAR_SCORE, sift=SIFT_SCORE, damage=DAMAGE[str(HVAR_PRED)])

				with open(ANN_ANNOVAR_AFTER, 'a') as f:
					f.write(RESULT)

		else:
			if GENE_DETAIL != "NA":
				if GENE_DETAIL.find('dist')!=-1:

					RESULT="{chr}_{pos}\t{ref}\t{alt}\t{var_cla}\t{dbsnp147}\t{iht}\t{func}\t{gene}\t{exonicfunc}\t-\t-\t{gene_detail}\t{spl}\t{trn}\t{hvar_score}\t{sift}\t{damage}\n".format(
						chr=CHR, pos=START, ref=REF, alt=ALT,
						var_cla=VAR_CLA, dbsnp147=DBSNP147, iht=IHT, func=FUNC,
						gene=GENE, exonicfunc=EXONICFUNC, gene_detail=GENE_DETAIL, spl=SPL,
						trn=TRN, hvar_score=HVAR_SCORE, sift=SIFT_SCORE, damage=DAMAGE[str(HVAR_PRED)])

					with open(ANN_ANNOVAR_AFTER, 'a') as f:
						f.write(RESULT)
				else:
					TRANS = GENE_DETAIL.split(':')[0]

					RESULT="{chr}_{pos}\t{ref}\t{alt}\t{var_cla}\t{dbsnp147}\t{iht}\t{func}\t{gene}\t{exonicfunc}\t{trans_id}\t-\t{gene_detail}\t{spl}\t{trn}\t{hvar_score}\t{sift}\t{damage}\n".format(
						chr=CHR, pos=START, ref=REF, alt=ALT,
						var_cla=VAR_CLA, dbsnp147=DBSNP147, iht=IHT, func=FUNC,
						gene=GENE, exonicfunc=EXONICFUNC, trans_id=TRANS, gene_detail=GENE_DETAIL,
						spl=SPL, trn=TRN, hvar_score=HVAR_SCORE, sift=SIFT_SCORE, damage=DAMAGE[str(HVAR_PRED)])

					with open(ANN_ANNOVAR_AFTER, 'a') as f:
						f.write(RESULT)
			else:
				RESULT="{chr}_{pos}\t{ref}\t{alt}\t{var_cla}\t{dbsnp147}\t{iht}\t{func}\t{gene}\t{exonicfunc}\t-\t-\t{aachange}\t{spl}\t{trn}\t{hvar_score}\t{sift}\t{damage}\n".format(
					chr=CHR, pos=START, ref=REF, alt=ALT,
					var_cla=VAR_CLA, dbsnp147=DBSNP147, iht=IHT, func=FUNC,
					gene=GENE, exonicfunc=EXONICFUNC, aachange=AACHANGE, spl=SPL,
					trn=TRN, hvar_score=HVAR_SCORE, sift=SIFT_SCORE, damage=DAMAGE[str(HVAR_PRED)])

				with open(ANN_ANNOVAR_AFTER, 'a') as f:
					f.write(RESULT)

def ClinvarMapping(ANNOVAR):

	# sqlite3
	# create table clinvar (chr_pos, ref, alt, mutation, aachange, rs, type, gene, class, pathogenic, disease, rcv, extra_id, other, pubmed)
	CLINVAR_ANN = '2.' + str(ANNOVAR).split(".")[2] + '.result'

	CONNECTION = sqlite3.connect(EXOME_ANN_DB)
	CUR = CONNECTION.cursor()

	with open(ANNOVAR, 'r') as ann:
		for line in ann.readlines():

			COL = line.strip().split('\t')

			CHR_POS = COL[0]
			REF = COL[1]
			ALT = COL[2]
			#TYPE = COL[3]
			#RS = COL[4]
			#RECESSIVE = COL[5]
			#REGION = COL[6]
			#GENE = COL[7]
			#FUNC = COL[8]
			TRANS = COL[9]
			#AACHANGE = COL[10]
			FULL_TRANS = COL[11]
			ETC = COL[12:]
			
			def PATHOGENIC (LINE):
				if LINE.find('benign') != -1:
					return 0
				elif LINE.find('Likely pathogenic')!= -1:
					return 1
				elif LINE.find('Pathogenic')!=-1:
					return 2
				else:
					return 0

			# search 
			if CUR.execute('select exists ( select * from clinvar where chr_pos = ? and ref = ? and alt = ?) ;', (CHR_POS, REF, ALT)).fetchone()[0] == 1:
				DB_COL = CUR.execute('select * from clinvar where chr_pos = ? and ref = ? and alt = ? ;', (CHR_POS, REF, ALT)).fetchone()

				DB_TRANS = DB_COL[3].split(':')[1]
				DB_CLASS = DB_COL[8]
				DB_DESIEASE = DB_COL[10]
				DB_EXTRA_ID = DB_COL[12]
				DB_PM =DB_COL[14]

				if TRANS == DB_TRANS:
					#RESULT = '\t'.join(COL.extend([PATHOGENIC(DB_CLASS), DB_DESIEASE,  DB_EXTRA_ID, DB_PM]))
					COL.extend([PATHOGENIC(DB_CLASS), DB_DESIEASE,  DB_EXTRA_ID, DB_PM])
					print("\t".join(COL))
					#print(COL)

				elif len(FULL_TRANS.split(','))!=1:
					for trans in FULL_TRANS.split(','):
						if trans.split(':')[1] == TRANS:
							AMINO_INFO = trans.split(':')[-1].split('.')[-1]
							COL[:11].extend(['p.' + AMINO_INFO[0]+AMINO_INFO[1:-1]+AMINO_INFO[-1], COL[12:], PATHOGENIC(DB_CLASS), DB_DESIEASE,  DB_EXTRA_ID, DB_PM])
							print("\t".join(COL))
							#print(COL)
							break
			else:
				COL.extend(['0','-','-','-'])
				print("\t".join(COL))
				#print(COL)

	CONNECTION.close()

def DriverMapping(CLINVAR):
	CONNECTION = sqlite3.connect(EXOME_ANN_DB)
	CUR = CONNECTION.cursor()

	with open(CLINVAR, 'r') as f:
		for line in f.readlines():
			COL=line.strip().split('\t')

			CHR = COL[0].split('_')[0]
			POS = COL[0].split('_')[-1]
			REF = COL[1]
			ALT = COL[2]

			if len(REF) == 1:
				if len(ALT) == 1:
					START = POS
					END = POS
				else:
					START = POS
					END = POS + len(ALT) -1
			else:
				if len(ALT) == 1:
					START = POS - len(REF)+1
					END = POS
				else:
					START = POS -len(REF) + 1
					END = POS + len(ALT) -1


			DB_COL = CUR.execute('select * from driver where chr = ? and pos >= ? and pos <= ? ;', (CHR, START, END)).fetchone()

			if DB_COL:
				if len(REF) == len(ALT):
					if REF == ALT:
						COL.extend(['2', DB_COL[-2], DB_COL[-1]])
					else:
						COL.extend(['1', DB_COL[-2], DB_COL[-1]])
				else:
					COL.extend(['1', DB_COL[-2], DB_COL[-1]])

			else: 
				COL.extend(['0', '-', '-'])

def snpEFF(ADI, DRIVERMAPPING):
	VCF = 
	ANN_VCF =
	ANN1_VCF= tempfile.NamedTemporaryFile(delete=False, bufsize=1)
	JAR_SNPEFF = '/tools/Bio-tools/snpEff/snpEff.jar'

	with open(VCF, 'w') as f:
		f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')

	with open(VCF, 'a') as f:
		with open(ADI, 'r') as adi:
			for line in adi:
				COL=line.strip().split()

				CHR=COL[1].split('_')[0]
				POS=COL[1].split('_')[1]
				REF=COL[4]
				ALT=COL[5]

				f.write('{chr}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.'.format(chr=CHR, pos=POS, ref=REF, alt=ALT))
	
	with open(ANN_VCF, 'w') as f:
		for line in subprocess.Popen(['java', '-jar', JAR_SNPEFF, '-v', 'hg19', VCF], stdout=subprocess.PIPE).stdout.readlines():
			f.write(line)

	with open(ANN1_VCF, 'w') as f:
		for line in subprocess.Popen(['sed', 's/||/|-|/g;s/||/|-|/g', ANN_VCF], stdout=subprocess.PIPE).stdout.readlines():
			f.write(line)

	with open(ANN_OUT, 'w') as ann_out:
		with open(ANN1_VCF, 'r') as f:
			for line in f.readlines():
				if line[0]!='#':
					COL=line.strip().split('\t')

					for INFO in COL[-1].split(','):

						ann_out.write('{chr}_{pos}\t{ref}\t{alt}\t{trans}\t{nnchange}\t{aachange}\t{gene}\t{impact}'.format(
							chr=COL[0], pos=COL[1], ref=COL[3], alt=COL[4],
							trans=INFO.split('|')[6], nnchange=INFO.split('|')[9], aachange=INFO.split("|")[10],
							gene=INFO.split('|')[3], impact=INFO.split('|')[2]))

	with open(SNPEFF_MAP, 'w') as snpeff_map:
		with open(ANN_OUT, 'r') as ann_out:
			for ann_out_line in ann_out:
				with open(DRIVERMAPPING, 'r') as drivermapping:
					for drivermapping_line in drivermapping:
						ANN_OUT_COL = ann_out_line.strip().split('\t')
						DRIVERMAPPING_COL = drivermapping_line.strip().split('\t')

						TRANS = DRIVERMAPPING_COL[11].split(':')[1]

						if DRIVERMAPPING_COL[11] != 'NA':
							if ANN_OUT_COL[0] ==  DRIVERMAPPING_COL[0] and ANN_OUT_COL[1] == DRIVERMAPPING_COL[1] and ANN_OUT_COL[2] == DRIVERMAPPING_COL[2] and ANN_OUT_COL[3].split('.')[0] == TRANS and DRIVERMAPPING_COL[10] == ANN_OUT_COL[5]:
								if DRIVERMAPPING_COL[16] == '-':
									DRIVERMAPPING_COL[:16].extend([ANN_OUT_COL[-1], DRIVERMAPPING_COL[17:]])
								else:
									DRIVERMAPPING_COL[:16].extend([DRIVERMAPPING_COL[16]+":"+ANN_OUT_COL[-1], DRIVERMAPPING_COL[17:]])







### MAIN ###
#def main ():

			

### ACTION ### 
if __name__ == "__main__":
	#argv = GetOpts()

	#main( )
	ADI=sys.argv[1]
	ClinvarMapping(ADI)

	sys.exit(0)
