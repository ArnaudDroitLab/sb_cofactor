import sys, glob, string, os

os.chdir("/home/chris/Bureau/sb_cofactor_hr/mESC")

def read_chrom_regions( filename ):

	chrom_regions = {}
	
	infile = open( filename, 'r' )
	
	lines = infile.readlines()
	
#	for line in lines[1:]:
	for line in lines:
#		line = string.split( line, '\t' )
		line = line.split('\t')
        
		try:
			chrom_regions[ line[0] ].append( [ int(line[1]), int(line[2]) ] )
		except KeyError:
			chrom_regions[ line[0] ] = [ [ int(line[1]), int(line[2]) ] ]
	
	infile.close
	
	for chrom in chrom_regions:
		chrom_regions[ chrom ].sort()

	return chrom_regions

def calculate_bp_length_chrom_regions( chrom_regions ):
	
	total_regions = 0
	total_bp = 0
	for chrom in chrom_regions:
		for region in chrom_regions[ chrom ]:
			total_regions += 1
			total_bp += region[1] - region[0]

	return total_regions, total_bp



def compare_two_chrom_regions( chrom_regions1, chrom_regions2 ):
	
	overlapped_bp = 0	
	for chrom1 in chrom_regions1:
		skip_count = 0
		for region1 in chrom_regions1[ chrom1 ]:
#			if chrom_regions2.has_key( chrom1 ):
			if chrom1 in chrom_regions2:
				for region2 in chrom_regions2[ chrom1 ][skip_count:]:
					if region1[0] <= region2[0] <= region1[1] or region2[0] <= region1[0] <= region2[1]:
						overlapped_bp += min( region1[1], region2[1] ) - max( region1[0], region2[0] )
					elif region2[0] > region1[1]:
						break
					elif region1[0] > region2[1]:
						skip_count += 1
	return overlapped_bp

def get_filename( filename_list ):
    file_list = []
    infile = open( filename_list, 'r' )
    lines = infile.readlines()
    for line in lines:
        path_bed = line.strip()
        file_list.append(path_bed)
    return(file_list)
        
genome_length = 2000000000

#file_list = glob.glob( 'ENRICHED_REGIONS*' )
# file_list = glob.glob( '*narrowPeak*' )
file_list = get_filename("/home/chris/Bureau/sb_cofactor_hr/mESC/scripts/framptongram/batch2/bed_file_path_v2_mergedreplicate_v2.txt")
list_chrom_regions = []
for filename in file_list:
	list_chrom_regions.append( read_chrom_regions( filename ) )

outfile = open( '/home/chris/Bureau/sb_cofactor_hr/mESC/scripts/framptongram/batch2/comparison_matrix_20190417_mergedreplicate.txt', 'w' )
# outfile.write( 'NAME' + '\t' + 'REGIONS' + '\t' + 'BP' + '\t' + string.join( file_list, '\t' ) + '\n' )
outfile.write( 'NAME' + '\t' + 'REGIONS' + '\t' + 'BP' + '\t' + '\t'.join(file_list) + '\n' )

for filename_counter, chrom_regions1 in enumerate( list_chrom_regions ):
	total_regions1, total_bp1 = calculate_bp_length_chrom_regions( chrom_regions1 )
#	print file_list[ filename_counter ]
	print(file_list[ filename_counter ])
	outfile.write( file_list[ filename_counter ] + '\t' + str( total_regions1 ) + '\t' + str( total_bp1 ) )
	for chrom_regions2 in list_chrom_regions:
		total_regions2, total_bp2 = calculate_bp_length_chrom_regions( chrom_regions2 )
		overlapped_bp = compare_two_chrom_regions( chrom_regions1, chrom_regions2 )
		score = ( float(overlapped_bp) / float(total_bp1)) / ( float(total_bp2) / float(genome_length) )
		mean1 = float( total_bp1 ) / genome_length
		mean2 = float( total_bp2 ) / genome_length
		stdev1 = pow( ( pow( 1 - mean1, 2) * float(total_bp1) + pow( mean1, 2) * float( genome_length - total_bp1 ) ) / float(genome_length), 0.5 )
		stdev2 = pow( ( pow( 1 - mean2, 2) * float(total_bp2) + pow( mean2, 2) * float( genome_length - total_bp2 ) ) / float(genome_length), 0.5 )
		correlation = ( float(overlapped_bp) - float(genome_length)* mean1 * mean2 ) / ( float(genome_length-1) * stdev1 * stdev2 )
		if correlation < 0:
			correlation = ( float(overlapped_bp+100) - float(genome_length)* mean1 * mean2 ) / ( float(genome_length-1) * stdev1 * stdev2 )
			if correlation < 0:
				min_correlation = ( float(0) - float(genome_length)* mean1 * mean2 ) / ( float(genome_length-1) * stdev1 * stdev2 )
				normalized_correlation = -pow( correlation / min_correlation, 0.85 )
			else:
				min_correlation = 0
				normalized_correlation = 0
#			print total_bp1, total_bp2, overlapped_bp, score, correlation, min_correlation, normalized_correlation
			print(total_bp1, total_bp2, overlapped_bp, score, correlation, min_correlation, normalized_correlation)
		elif correlation > 0:
			max_correlation = ( float( min( total_bp1, total_bp2 ) ) - float(genome_length)* mean1 * mean2 ) / ( float(genome_length-1) * stdev1 * stdev2 )
			normalized_correlation = pow( correlation / max_correlation, 0.85 )
#			print total_bp1, total_bp2, overlapped_bp, score, correlation, max_correlation, normalized_correlation
			print(total_bp1, total_bp2, overlapped_bp, score, correlation, max_correlation, normalized_correlation)


		outfile.write( '\t' + str( normalized_correlation ) )
	outfile.write( '\n' )\


outfile.close()
