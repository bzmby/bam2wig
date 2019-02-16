#Author: Behzad Moumbeini
#how to run: python bam_to_wig.py -o outfile infile.bam

import os
import HTSeq
import argparse


def empty_array_from_file(bam_file,stranded=True, typecode="i"):
	cov_array=HTSeq.GenomicArray("auto", stranded = stranded, typecode = typecode)	
	myheader = HTSeq.BAM_Reader(bam_file).get_header_dict()
	for entry in myheader['SQ']:
		cov_array.add_chrom(entry['SN'],entry['LN'])
	return cov_array
    



#Read BAM files and load the alignments on a coverage track
def load_array (Gen_array,infile):
	print "loading data from %s" %(infile)
	for algn in HTSeq.BAM_Reader(infile):
		if algn.aligned:
			if algn.iv.chrom in Gen_array.chrom_vectors.keys():
				for algn_cig in algn.cigar:
					if algn_cig.type == 'M':
						Gen_array[algn_cig.ref_iv] += 1
					else:
						continue


#Write a track in wiggle format from an unstranded GenomicArray and normalize by a factor
def write_wig_track(cov_array,outprefix,norm_factor=1):
	strands = ['+', '-'] if cov_array.stranded else ["."]
	for strand in strands:
		outf = open("%s_track_%s_.wig" %(outprefix,strand), 'w')
		print "Writing file %s_track_%s_.wig..." %(outprefix,strand)
		if strand == '+':
			mycol = '255,153,51'
		elif strand == '-':
			mycol = '51,102,153'
		else:
			mycol = '0,0,0'
		outf.write("track type=wiggle_0 alwaysZero=on name=%s_%s visibility=full maxHeightPixels=30 graphType=bar color=%s\n" %(os.path.basename(outprefix),strand, mycol))
		for vect in [i[strand] for i in cov_array.chrom_vectors.values()]:
			tmp_vals=list()
			start = ''
			for iv,val in vect.steps():
				if val == 0:
					if (not tmp_vals):
							continue
					elif iv.length <= 10:
							tmp_vals.extend([val]*iv.length)
					elif iv.length > 10:
						outf.write("fixedStep chrom=%s start=%d step=1\n" %(iv.chrom,start))
						for i in tmp_vals:
							outf.write("%.6f\n" %(i*norm_factor))
						tmp_vals=list()
						start=''
				elif val != 0 and (not tmp_vals):
					start = iv.start + 1
					tmp_vals.extend([val]*iv.length)
				elif val !=0 and tmp_vals:
					tmp_vals.extend([val]*iv.length)
			if tmp_vals:
				outf.write("fixedStep chrom=%s start=%d step=1\n" %(iv.chrom,start))
				for i in tmp_vals:
					outf.write("%.6f\n" %(i*norm_factor))
		outf.close()
		print "File  %s_track_%s_.wig sucessfully created" %(outprefix, strand)



def main():
	parser=argparse.ArgumentParser(description='Makes a coverage track from unpaired stranded bam files')
	parser.add_argument('infile',help='BAM file')
	parser.add_argument('-o','--outf', help='Output filename', default='coverage_track')
	args=parser.parse_args()

	
	#Output filename
	outfile=args.outf
	if os.path.dirname(outfile) and (not os.path.exists(os.path.dirname(outfile))):
		os.makedirs(os.path.dirname(outfile))
	
	#Do the job
	cvg=empty_array_from_file(args.infile) #Create an empty genomic array from file
	load_array(cvg,args.infile) #Loads the array with data
	write_wig_track(cvg, outfile) #Write the result in wig format

if __name__ == '__main__':
    main()
