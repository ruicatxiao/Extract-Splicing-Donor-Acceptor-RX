import subprocess

# Load AGAT and convert GFF file
subprocess.run(['ml', 'AGAT/0.8.1'])
subprocess.run([
    'agat_convert_sp_gxf2gxf.pl', '-v', '3',
    '-g', 'WebApollo_SN3e1_230701_raw.gff3',
    '-o', 'WebApollo_SN3e1_230701_cleaned.gff'
])

# Add introns to GFF file
subprocess.run([
    'agat_sp_add_introns.pl',
    '--gff', 'WebApollo_SN3e1_230701_cleaned.gff',
    '--out', 'mrna_apollo_intron_added.gff'
])

# Extract intron lines
with open('mrna_apollo_intron_added.gff', 'r') as infile, \
     open('WebApollo_SN3e1_230701_intron_only.gff', 'w') as outfile:
    for line in infile:
        if line.split()[2] == "intron":
            outfile.write(line)

# Process intron data and generate intron.bed
intron_data = []
with open('WebApollo_SN3e1_230701_intron_only.gff', 'r') as infile:
    for line in infile:
        fields = line.strip().split('\t')
        attributes = dict(attr.split('=') for attr in fields[8].split(';'))
        parent_id = attributes['Parent']
        id = attributes['ID']
        intron_data.append((fields[0], fields[3], fields[4], parent_id, id, fields[6]))

with open('intron.bed', 'w') as outfile:
    for chrom, start, end, parent, id, strand in intron_data:
        start = int(start) - 10
        end = int(end) + 10
        outfile.write(f"{chrom}\t{start}\t{end}\t{parent}_{id}\t.\t{strand}\n")

# Process donor and acceptor data and generate donor and acceptor bed files
intron_donor_data = []
intron_acceptor_data = []

with open('intron.bed', 'r') as infile:
    for line in infile:
        fields = line.strip().split('\t')
        chrom, start, end, identifier, _, strand = fields
        if strand == "+":
            intron_donor_data.append((chrom, str(int(start) - 1), str(int(start) + 21), f"DONOR_{identifier}", ".", strand))
        elif strand == "-":
            intron_donor_data.append((chrom, str(int(end) - 22), end, f"DONOR_{identifier}", ".", strand))
            intron_acceptor_data.append((chrom, str(int(start) - 1), str(int(start) + 21), f"ACCEPTOR_{identifier}", ".", strand))
        intron_acceptor_data.append((chrom, str(int(end) - 22), end, f"ACCEPTOR_{identifier}", ".", strand))

# Sort donor and acceptor data and write to files
intron_donor_data.sort(key=lambda x: (x[0], int(x[1])))
intron_acceptor_data.sort(key=lambda x: (x[0], int(x[1])))

with open('intron_donor.bed', 'w') as outfile:
    for data in intron_donor_data:
        outfile.write('\t'.join(data) + '\n')

with open('intron_acceptor.bed', 'w') as outfile:
    for data in intron_acceptor_data:
        outfile.write('\t'.join(data) + '\n')

# Use bedtools to extract sequences
subprocess.run([
    'bedtools', 'getfasta', '-s',
    '-fi', 'sn3_hic_genome.fa',
    '-fo', 'intron_donor.fa',
    '-bed', 'intron_donor.bed',
    '-nameOnly'
])

subprocess.run([
    'bedtools', 'getfasta', '-s',
    '-fi', 'sn3_hic_genome.fa',
    '-fo', 'intron_acceptor.fa',
    '-bed', 'intron_acceptor.bed',
    '-nameOnly'
])
