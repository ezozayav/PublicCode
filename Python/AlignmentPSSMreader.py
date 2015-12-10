#Sep 22 2014
#dr.mark.schultz@gmail.com
#github: schultzm
#This file will return an alignment pssm, with SNPs called, from an input PSSM file.  
#this does not work properly.  will not write the binary pssm
#First, import the modules
import argparse

# sys.exit(0)
parser = argparse.ArgumentParser(description = "This file will read and format a PSSM file for plotting in R.")
parser.add_argument('-i', '--input', help = 'Input PSSM file to scan.', required=True)
parser.add_argument('-g', '--gaps_as_fifth_state', help = 'Count gaps as a fifth state?', default = 'No', required = False)

args = parser.parse_args()

def PSSM_input(input):
	with open(input, 'r') as input_handle:
		next(input_handle)
		if args.gaps_as_fifth_state[0].lower() == 'n':
			with open(input+'PSSM_variantSitesCalled_gapsNotCounted.txt', 'w') as output_handle:
				output_handle.write('Ref\tGap\tA\tC\tG\tT\tSNP\tPyIndex\n')
				pos = 0
				for line in input_handle:
					pssm_list = line.rstrip('\n').split('\t')
					binary_pssm = [1  if float(i) > 0 else 0 for i in pssm_list[3:]]
					print pssm_list
					print binary_pssm
					if sum(binary_pssm) > 2:
						pssm_list.append('variant')
					else:
						pssm_list.append('invariant')
					pssm_list.append(pos)
					output_handle.write(str(('\t').join(str(x) for x in pssm_list))+'\n')
					pos += 1
		else:
			#gaps counted as a fifth state here
			with open(input+'PSSM_variantSitesCalled_gapsFifthState.txt', 'w') as output_handle:
				output_handle.write('Ref\tGap\tA\tC\tG\tT\tSNP\tPyIndex\n')
				pos = 0
				for line in input_handle:
					pssm_list = line.rstrip('\n').split('\t')
					binary_pssm = [1  if float(i) > 0 else 0 for i in pssm_list[2:]]
#					print pssm_list
#					print binary_pssm
					if sum(binary_pssm) > 2:
						pssm_list.append('variant')
					else:
						pssm_list.append('invariant')
					pssm_list.append(pos)
					output_handle.write(str(('\t').join(str(x) for x in pssm_list))+'\n')
					pos += 1

PSSM_input(args.input)
