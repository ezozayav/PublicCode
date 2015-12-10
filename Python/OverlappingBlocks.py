#Sep 23 2014
#dr.mark.schultz@gmail.com
#github: schultzm
#This file will pull out the merged-overlapping alignment coordinates in a table of coordinates with start column and stop column


#First, import the modules
import argparse

parser = argparse.ArgumentParser(description = "This file will pull out the merged-overlapping alignment coordinates in a table of coordinates with start column and stop column")
parser.add_argument('-i', '--input_block_coordinates', help = 'Input csv which contains the block coordinates.', required=True)
parser.add_argument('-p', '--prefix', help = 'Type something meaningful here to allow tracing of the file trail', required=True)

args = parser.parse_args()

#pull the alignments from the cluster numbered 'cluster_number' to get coordinates of the supermatrices 
def filter(input, prefix):
	with open(input, 'r') as input_handle:
		next(input_handle)
		alignments=[]
		for line in input_handle:
			x = line.rstrip('\n').split(",")
			alignments.append(x)
		#alignments is the list containing coordinates for downstream extraction
		blocks = []
		ranges = []
		merged_blocks = []
		for i in alignments:
			#append this mini list "i" to the large blocks list
			i = map(int, i)
			blocks.append(i)
		blocks.sort()
# 		for i in blocks:
# 			print type(i)
# 			for j in i:
# 				print type(j)
# 		print type(blocks)
		with open(prefix+'_BlocksSigAlignment_coordinates.txt', 'w') as output_handle:
			output_handle.write('start, stop\n')
			for group in yield_data(blocks):
				towrite = str(group[0])+', '+str(group[1])+'\n'
				output_handle.write(towrite)
# 		print blocks

#this function will merge the blocks		
def yield_data(data):
	start = None
	end = None
	for entry in sorted(data):
		if start is None:
			start = entry[0]
		if end is None or entry[0] <= end:
			end = entry[1]
		else:
			yield [start, end]
			start, end = entry
	yield [start, end]

filter(args.input_block_coordinates, args.prefix)
print 'Done'
