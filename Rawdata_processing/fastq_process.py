from Bio import SeqIO
from gzip import open as gzopen

input_filename = "PATHTOREAD2/R2.fastq.gz"
output_filename = "sample_R2_processed.fastq"

output_handle = open(output_filename, "w")

with gzopen(input_filename, "rt") as handle:
    for record in SeqIO.parse(handle, "fastq"):
      cut_record = record[32:40] + record[70:78] + record[22:32]  # BC2 + BC1 + UMI
      SeqIO.write(cut_record, output_handle, "fastq")

output_handle.close()