#! /usr/bin/python3 -u

# ##############################################################################
# Performs parsing of the cleaned fastq files to identify the barcodes.
# Written and designed by Jose Rojas Echenique
# Modified by Alex N Nguyen Ba and Katherine Lawrence
# rojasechenique@fas.harvard.edu and nnguyenba@fas.harvard.edu
# 2017
#
# Version 1
#
# To run, use clean script on the Fastq file.
# Then run this script on the cleaned fastq file.
#
# LICENCE 
#
# ##############################################################################

import sys
from collections import OrderedDict
#sys.path.append('/n/desai_lab/users/klawrence/BBQ/pipeline')

import regex

complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def rc(seq):
	return "".join(complements.get(base, base) for base in reversed(seq))

#def make_corrector(options):
#	checkers = [regex.compile("("+o+"){e<=2}") for o in options]
#	def corrector(match):
#		for (i,c) in enumerate(checkers):
#			if c.match(match):
#				return i+1
#	return corrector

#Alex's new corrector
def make_corrector(options, num_mismatch=0):
	checkers = [regex.compile("("+o+"){e<=" + str(num_mismatch) + "}") for o in options]
	def corrector(match):
		current_error_count = 10000
		current_best_i = ""
		for (i,c) in enumerate(checkers):
			m = c.fullmatch(match)
			if m:
				error_count = m.fuzzy_counts[0] + m.fuzzy_counts[1] + m.fuzzy_counts[2]
				if(error_count == 0):
					return i+1
				elif (error_count < current_error_count and error_count <= num_mismatch):
					current_best_i = i+1
					current_error_count = error_count
		if (current_best_i != ""):
			return current_best_i
		else:
			return 0
	return corrector


molecular_barcode = ".{8}"
lineage_barcode_BY = "(.{3})(?e)(?:GC){e<=1}(.{5})(?:AT){e<=1}(.{5})(?:AT){e<=1}(.{3})"
lineage_barcode_RM = "(.{3})(?e)(?:AT){e<=1}(.{5})(?:GC){e<=1}(.{5})(?:AT){e<=1}(.{3})"
lineage_barcode_BYrc = "(.{3})(?e)(?:AT){e<=1}(.{5})(?:AT){e<=1}(.{5})(?:GC){e<=1}(.{3})"
lineage_barcode_RMrc = "(.{3})(?e)(?:AT){e<=1}(.{5})(?:GC){e<=1}(.{5})(?:AT){e<=1}(.{3})"

#lineage_barcode = "(.{18,22})"

KanP1 = "(?e)(?:GCAATTAACCCTCACTAAAGG)"
HygP1 = "(?e)(?:CCGCTGAGCAATAACTAGC)"
HygP1rc = "(?e)(?:GCTAGTTATTGCTCAGCGG)"

P5modindex = [
	"TACGAG",
	"CCAATG",
	"GACCTT",
	"TGCAGT",
	"TCTTGT",
	"AAGTGT",
	"CCTTCG",
	"GCAAGA",
	"AATACC",
	"AGCGTT",
	"CTGATT",
	"AAGGCG",
	"TGAACG",
	"TTGCTG",
	"TATCCT",
	"GGTACT",
	"TTCTCC",
	"TTGGAC",
	"TTACGC",
	"TCCGTA",
	"CTTCAC",
	"CTTAGG",
	"CTCTGA",
	"CGTGAA",
	"CGACTA",
	"CAGTTC",
	"GTTCGA",
	"GTCAAC",
	"GTAGTC",
	"GCTCAT",
	"GCCTAA",
	"GAGCAA",
	"GAATGG",
	"ACTGTG",
	"ACGTCA",
	"TTATAA"]

correct_P5mod = make_corrector(P5modindex,1)

P5indexDict = {idx:(i+1) for i,idx in enumerate(P5modindex)}

bbq_col_index = [
	"TTGATCG",
	"GCATTGAT",
	"AAGGCGTCT",
	"CATTCTCGAG",
	"CCTTAGTACTA",
	"AGTGGAATACCT",
	"GGTATCA",
	"GATGTCCT",
	"CTCATCCAG",
	"GACGGAACTC",
	"ATCGCAGGCAT",
	"TCACGACTAGTA"]


correct_bbq_col = make_corrector(bbq_col_index,1)

bbq_row_index = [
	"TACCTGA",
	"AGAACCAT",
	"ACCGTTGAC",
	"GCAAGCTGTT",
	"TCGGTGGTACG",
	"GTTCAACCGATT"]

correct_bbq_row = make_corrector(bbq_row_index,1)

#ilt_row_index_fwd = [rc(x) for x in ilt_row_index]
#correct_ilt_row_fwd = make_corrector(ilt_row_index_fwd[0:7])

P7modindex = [
	"AGCTCC",
	"TTATGG",
	"CGTTAC",
	"TTGTCC",
	"ACAGTC",
	"CGGTGA",
	"TTGGAA",
	"CCTGGA",
	"ATCAGG",
	"ACAACT",
	"GCACAT",
	"CGTCTT",
	"TGTGGC",
	"GGATTG",
	"TCGGTT",
	"ATGCAC",
	"AGTTGG",
	"TCTCAG",
	"TAGATC",
	"CACTCT",
	"AACCAA",
	"GTTGTT",
	"CTAGGT",
	"GAGTAG",
	"AAGGCT",
	"AAGTTA",
	"AATATT",
	"AACGTG",
	"AGTGAA",
	"ATCTAT",
	"GAAGCG",
	"GAATGA",
	"GGTCCA",
	"TTCCTC",
	"GTCATA",
	"GTCGAC"]

correct_P7mod = make_corrector(P7modindex,1)

P7indexDict = {idx:(i+1) for i,idx in enumerate(P7modindex)}

#new lox site
Lox = "(?e)(?:ATAACTTCGTATAATGTATGCTATACGAAGTTAT)"
Loxrc = "(?e)(?:ATAACTTCGTATAGCATACATTATACGAAGTTAT)"

#old lox sites
#LoxP2R = "(?e)(?:GAAAGGGGCTATAGCATACATTATAGCCCCTTTC)"
#LoxP2Rrc = "(?e)(?:GAAAGGGGCTATAATGTATGCTATAGCCCCTTTC)"
#Lox5171L = "(?e)(?:GAAAGGGGCTATAATGTGTACTATACGAAGTTAT)"
#Lox5171Lrc = "(?e)(?:ATAACTTCGTATAGTACACATTATAGCCCCTTTC)"

#Lox5171RL = "(?e)(?:GAAAGGGGCTATAATGTGTACTATAGCCCCTTTC)"
#Lox5171RLrc = "(?e)(?:GAAAGGGGCTATAGTACACATTATAGCCCCTTTC)"

#pTEFRrc = "(?e)(?:GGGGACGAGGCAAGCT)"
#LoxP2Lrc = "(?e)(?:ATAACTTCGTATAATGTATGCTATAGCCCCTTTC)"
#LoxPLRrc = "(?e)(?:GAAAGGGGCTATAATGTATGCTATAGCCCCTTTC)"
#pTEFmidRrc = "(?e)(?:GGAGCCGTAATTTTTGCTTC)"

def OR(xs):
	return "(" + "|".join(["(?:"+x+")" for x in xs]) + ")"


BC_R1 = (
	    "(" + molecular_barcode + ")" +
	    OR(bbq_col_index) + "{e<=1}" +
	    KanP1 + "{e<=4}" +
		 '(?e)(?:AGCG){e<=1}' + 
	    lineage_barcode_BY +
		 '(?e)(?:AGTA){e<=1}' )
		 

BC_R2 = (
		"(" + molecular_barcode + ")" +
		OR(bbq_row_index) + "{e<=1}" +
		HygP1rc + "{e<=4}" +
		'(?e)(?:TACT){e<=1}' +
		lineage_barcode_RMrc + 
		'(?e)(?:CGCT){e<=1}')
		 
#S5 = (OR(S5index))
#N7 = (OR(N7index))	

BC_R1re = regex.compile(BC_R1)
BC_R2re = regex.compile(BC_R2)
#S5re = regex.compile(S5)
#N7re = regex.compile(N7)


#KanP1re = regex.compile((OR(ilt_col_index[0:12]) + "{e<=2}" + KanP1 + "{e<=4}"))
#HygP1re = regex.compile((HygP1 + "{e<=4}" + OR(ilt_row_index_fwd[0:7]) + "{e<=2}"))
#LoxP2Rre = regex.compile((LoxP2R + "{e<=4}"))
#BCre = regex.compile(('TACT'+ "{e<=1}" + lineage_barcode +'CGCT'+ "{e<=1}"))
#N7_R1re = regex.compile((OR(N7index)))
#N7rc_R1re = regex.compile((OR(N7rc)))
#S5_R1re = regex.compile((OR(S5index)))
#S5rc_R1re = regex.compile((OR(S5rc)))


nreads = 0
n_match = 0
r1_discarded = 0
r2_discarded = 0
for line in sys.stdin:
	nreads += 1
	row=OrderedDict([
		    (x,"") for x in ["illumina_index_P5","illumina_index_P7","col","row","barcode_1","barcode_2","mb1","mb2"]
	    ])
	
	splitline = line.split("	")
	if len(splitline) < 4:
		continue

	#R2 = splitline[0]
	#R3 = splitline[1]
	R1 = splitline[2]
	R2 = splitline[3]
 
	m1 = BC_R1re.match(R1)
	if m1:
		m2 = BC_R2re.match(R2)
		if m2:
			n_match += 1
			row["illumina_index_P7"] = correct_P7mod(splitline[0])
			row["illumina_index_P5"] = correct_P5mod(splitline[1])
			row["mb1"] = m1.groups()[0]
			row["col"] = correct_bbq_col(m1.groups()[1])
			row["barcode_1"] = m1.groups()[2]+m1.groups()[3]+m1.groups()[4]+m1.groups()[5]
			row["barcode_2"] = rc(m2.groups()[2]+m2.groups()[3]+m2.groups()[4]+m2.groups()[5])
			row["row"] = correct_bbq_row(m2.groups()[1])
			row["mb2"] = m2.groups()[0]
			
			print(*row.values(), sep="	", file=sys.stdout)	

		else:
			r2_discarded += 1
			continue
	else:
		r1_discarded += 1
		continue
 
	
 
print("Discarded " + str(r1_discarded) + " + " + str(r2_discarded) + " out of " + str(nreads) + " reads.", file=sys.stderr)
print("Parsed "+str(n_match/nreads*100)+"% of reads",file=sys.stderr)


