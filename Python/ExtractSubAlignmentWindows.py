#May 14 2014
#dr.mark.schultz@gmail.com
#github: schultzm
#This file will split an alignment into overlapping sub alignments and output each sub-alignment to a file, and will output the alignment remainder too.


#First, import the modules
import argparse
import sys
from Bio import AlignIO


#set up the arguments parser to deal with the command line input
parser = argparse.ArgumentParser(description = "This script extracts sub alignments from a DNA alignment (might work for protein alignments; however, this hasn't been tested).  Input can be any alignment format. Output can only be fasta or phylip (to do: fix script to output any format).  All sub-alignments will be of length '-w', until the end of the alignment is reached.  The last sub-alignment will be the remainder.  Sub-alignment extractions will be made at each slide position, moving 5-prime to 3-prime.  The size of the overlap between sub-alignments equals 'window minus slide'.  The slide interval must be less than the window size, and the window or the slide must be smaller than the alignment length.  Note: for SNP alignments where the coordinates of the SNPs in the genome are no longer known, the coordinates of the window in the output file name will be meaningless when making comparisons to the original whole-genome sequence.  In this situation, one would need to compare the window coordinates of the output file to the original SNP file that contains the genome coordinates.  For help, email the author at dr.mark.schultz@gmail.com.")
parser.add_argument('-i', '--input', help = 'Input alignment file to pull the sub-alignments from.', required=True)
parser.add_argument('-f', '--informat', help = 'Informat of input file.', required = True)
parser.add_argument('-o', '--outformat', help = 'Outformat of output file (fasta or phylip only). NB: if specifying Phylip outformat, input-file sequence names must be <= 10 characters.', required = True)
parser.add_argument('-w', '--windowsize', type = int, help = 'Size of window (sites) output to each file.', required = True)
parser.add_argument('-s', '--slide', type = int, help = 'Slide this number of sites to direction 3-prime on subsequent sub-alignment extractions.', required = True)
args = parser.parse_args()

#list the command line arguments as input by the user
print '\nFile from which sub-alignments will be extracted: %s' % args.input,'\n'
print 'File informat: %s' % args.informat
print 'File outformat: %s' % args.outformat
print 'Window size: %s' % args.windowsize, 'contiguous sites'
if args.windowsize == 0:
	print "\nError: Window must be greater than zero!  At this setting, zero sites would be sampled. \nProgram exiting.\n"
	sys.exit(0)
print 'Slide: %s' % args.slide, 'sites to 3-prime'
print 'Slide interval as a percentage of window size: %s ' % (100*(float(args.slide)/float(args.windowsize))), 'percent'
print 'Note: the overlap between a sub-alignment and its previous sub-alignment will be %s' % (args.windowsize-args.slide), 'sites\n'


def split_alignment(input, informat, outformat, windowsize, slide):
	alignment = AlignIO.read(input, informat)
	length = alignment.get_alignment_length()
	print 'Alignment length:',length,'sites.\n'
	windowsize=int(windowsize)
	#error checks
	if windowsize-length > 0:
		print 'Error: Window size',windowsize,'is', windowsize-length, 'site(s) larger than input alignment length! \nProgram exiting.\n'
	elif slide-length > 0:
		print 'Error: Slide',slide,'is', slide-length, 'sites(s) larger than input alignment length! \nThis would result in an unwanted move to 5-prime.  \nRe-run whilst ensuring that your slide is less than or equal to your alignment length.  \nProgram exiting.\n'
	elif slide - windowsize > 0:
		print 'Error: Slide is greater than window size!  \nThis would result in sections of the alignment being skipped at subsequent window slides.  \nRe-run whilst ensuring that your slide is less than or equal to your window.  \nProgram exiting.\n'
	elif slide == 0:
		print "Error: Slide must be greater than zero!  At this setting, the same window would be sampled an unending number of times. \nProgram exiting.\n"  
	#run the function
	else:
		print 'Subsampling input alignment:'
		origin=0
		#subsample alignment in window sized blocks, sliding along each time by the slide value
		while origin+windowsize < length:
			with open(str(input.replace('.'+str(informat),''))+'_outputCoords_site'+str(origin+1)+'to'+str(origin+windowsize)+'_'+str((origin+windowsize)-origin)+'sitesInclusive.'+str(outformat), 'w') as handle:
				#check section 6.3.1 'Slicing alignments' in http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec84 to make sense of this next line
				AlignIO.write(alignment[:, origin:origin+windowsize], handle, outformat)
				handle.close()
			print "writing sites",origin+1,":",origin+windowsize,"to file '"+str(input.replace("."+str(informat),""))+"_outputCoords_site"+str(origin+1)+"to"+str(origin+windowsize)+"_"+str((origin+windowsize)-origin)+"sitesInclusive."+str(outformat)+"'"
			#error check will terminate the script if slide rounds to zero
			if int(float(windowsize)*(float(slide)/float(windowsize))) < 1:
				print 'slide rounded:',int(float(windowsize)*(float(slide)/float(windowsize)))
				print '\nError: Strangely, this combination of window and slide has resulted in the slide value being rounded to zero!  \nAt this setting combination, the same window would be sampled an unending number of times.  \nChange the window or slide setting and re-run. \nProgram exiting.\n'
				sys.exit(0)
			origin+=int(float(windowsize)*(float(slide)/float(windowsize)))
		#grab the remainder of the alignment whose length will be less than the window size
		else:
			with open(str(input.replace('.'+str(informat),''))+'_outputCoords_site'+str(origin+1)+'to'+str(length)+'_'+str(length-origin)+'sitesInclusive.'+str(outformat), 'w') as handle:
				#check section 6.3.1 'Slicing alignments' in http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec84 to make sense of this next line
				AlignIO.write(alignment[:, origin:length], handle, outformat)
				handle.close()
			print "writing sites",origin+1,":",length,"to file '"+str(input.replace('.'+str(informat),''))+'_outputCoords_site'+str(origin+1)+'to'+str(length)+'_'+str(length-origin)+'sitesInclusive.'+str(outformat)+"'"
			origin+=windowsize*20000000
			print "\nFinished.\n"

#on the input file, execute the split_alignment function
split_alignment(args.input, args.informat, args.outformat, args.windowsize, args.slide)
