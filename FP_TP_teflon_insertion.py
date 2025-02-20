import argparse

# Function to parse a BED file
def parse_bed(bed_file):
    bed_data = []
    with open(bed_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom, start, end, name1, name2, strand = parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4], parts[5]
            bed_data.append((chrom, start, end, name1, name2, strand))
    return bed_data

# Function to parse a VCF file
def parse_vcf(vcf_file):
    vcf_data = set()
    with open(vcf_file, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                parts = line.strip().split('\t')
                chrom, pos, name, ref, alt, quality,info = parts[0], int(parts[1]), parts[2], parts[3], parts[4], parts[6], parts[7]
                vcf_data.add((chrom, pos, name, ref, alt, quality, info))
    return vcf_data


parser = argparse.ArgumentParser(description="Compare the insertions detected by TEFLoN to a vcf from variant calling analysis on the same sample and on the reference annotation file")
parser.add_argument('-nbed', dest="nbed_file", required=True, help="Input BED file with non-reference insertions")
parser.add_argument('-rbed', dest="rbed_file", required=True, help="Input BED file with reference insertions")
parser.add_argument('-vcf', dest="vcf_file",required=True, help="Input VCF file")
parser.add_argument('-ref', dest="ref_file",required=True, help="Input reference BED file")
parser.add_argument("-prefix", required=True, help="Sample prefix for output file names")

args = parser.parse_args()

######## Classification of the insertions not present in the reference genome ########
# Load the BED and VCF data
newbed_data = parse_bed(args.nbed_file)
vcf_data = parse_vcf(args.vcf_file)

# Output files
new_true_positives_bed = open(f"{args.prefix}_new_true_positives.bed", "w")
new_false_positives_bed = open(f"{args.prefix}_new_false_positives.bed", "w")
new_false_negatives_vcf = open(f"{args.prefix}_new_false_negatives.vcf", "w")

# Compare non-reference insertions and classify entries
for chrom, start, end, name, family, strand in newbed_data:
    found = False
    for vcf_chrom, vcf_pos, vcf_name, ref, alt, quality, info in vcf_data:
        if chrom == vcf_chrom and (vcf_pos >= start - 20 and vcf_pos <= end):
            new_true_positives_bed.write(f"{chrom}\t{start}\t{end}\t{vcf_name}\t{family}\t{strand}\n")
            found = True
            break

    if not found:
        new_false_positives_bed.write(f"{chrom}\t{start}\t{end}\t{name}\t{family}\t{strand}\n")

for vcf_chrom, vcf_pos, vcf_name, ref, alt, quality, info in vcf_data:
    found = False
    for chrom, start, end, name, family, strand in newbed_data:
        if chrom == vcf_chrom and (vcf_pos >= start - 20 and vcf_pos <= end + 20):
            found = True
            break

    if not found:
        new_false_negatives_vcf.write(f"{vcf_chrom}\t{vcf_pos}\t{vcf_name}\t{ref}\t{alt}\t.\t{quality}\t{info}\n")

# Close output files
new_true_positives_bed.close()
new_false_positives_bed.close()
new_false_negatives_vcf.close()

######## Classification of the insertions also prensent in the reference genome ########
# Load the BED and ref data
refbed_data = parse_bed(args.rbed_file)
ref_data = parse_bed(args.ref_file)

# Output files
ref_true_positives_bed = open(f"{args.prefix}_ref_true_positives.bed", "w")
ref_false_positives_bed = open(f"{args.prefix}_ref_false_positives.bed", "w")

# Compare ref insertions and classify entries
for chrom, start, end, name, family, strand in refbed_data:
    found = False
    for ref_chrom, ref_start, ref_end, ref_name, dot, ref_strand in ref_data:
        if chrom == ref_chrom and name == ref_name and ((ref_end >= start >= ref_start - 20) or (end >= ref_start >= start)):
            ref_true_positives_bed.write(f"{chrom}\t{start}\t{end}\t{name}\t{family}\t{strand}\n")
            found = True
            break

    if not found:
        ref_false_positives_bed.write(f"{chrom}\t{start}\t{end}\t{name}\t{family}\t{strand}\n")

ref_true_positives_bed.close()
ref_false_positives_bed.close()

