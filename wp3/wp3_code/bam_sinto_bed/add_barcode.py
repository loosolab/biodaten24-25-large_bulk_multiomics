import os
import sys
import pysam
# Adapted from https://www.biostars.org/p/450160/
def add_barcode_to_bam(input_bam, barcode, outfile):
    output_bam = outfile # input_bam.replace(".bam", ".cb.bam")
    # The output bam needs to be in OUTPUT_BAM_BC
    # Chane path accordingly
    
    infile = pysam.AlignmentFile(input_bam, "rb")
    outfile = pysam.AlignmentFile(output_bam, "wb", header=infile.header)
    
    for read in infile:
        read.set_tag("CB", barcode, "Z", replace=True)
        outfile.write(read)
        
    infile.close()
    outfile.close()

if __name__ == "__main__":   
    input_bam = sys.argv[1]
    barcode = sys.argv[2]
    outfile = sys.argv[3]
    
    add_barcode_to_bam(input_bam, barcode, outfile)
