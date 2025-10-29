import gzip

# Define the input and output file names
input_file = 'ensp_pfam_v115.domtblout.txt.gz'
output_file = 'results/isoform_pfam_domains.txt'

# The specific header line to be kept from the original file
header_to_keep = '# target name'

# A flag to ensure the header is written only once
header_written = False

# Open the gzipped input file for reading in text mode, and the output file for writing
with gzip.open(input_file, 'rt') as f_in, open(output_file, 'w') as f_out:
    # Iterate over each line in the input file
    for line in f_in:
        # Check if the line is the header we want to keep
        if line.strip().startswith(header_to_keep):
            # If the header hasn't been written yet, write it to the output file
            if not header_written:
                f_out.write(line)
                header_written = True
        # Check if the line is a data line (i.e., it doesn't start with '#')
        elif not line.strip().startswith('#') and line.strip():
            # Write the data line to the output file
            f_out.write(line)

print(f"Cleaned data has been saved to {output_file}")
