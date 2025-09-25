import os
import pysam
from intervaltree import IntervalTree
import glob
import multiprocessing

def read_variants(vcf_path):
    variants = []
    with pysam.VariantFile(vcf_path) as vcf_in:
        for record in vcf_in:
            variants.append({
                'chrom': record.chrom,
                'start': record.start,
                'end': record.info.get('END', record.stop),
                'svtype': record.info.get('SVTYPE', None),
                'qual': record.qual
            })
    return variants

def build_interval_trees_by_chrom(variant_list):
    trees = {}
    for var in variant_list:
        chrom = var['chrom']
        start = var['start']
        end   = var['end']
        if chrom not in trees:
            trees[chrom] = IntervalTree()
        trees[chrom][start:end] = (start, end, var.get('svtype'))
    return trees

def find_overlaps_in_reference(de_novos, ref_trees, overlap_threshold=70): 
    found = []
    not_found = []
    for var in de_novos:
        chrom = var['chrom']
        if chrom not in ref_trees:
            not_found.append(var)
            continue
        intervals = ref_trees[chrom].overlap(var['start'], var['end'])
        is_found = False
        for interval in intervals:
            ref_start, ref_end, ref_type = interval.data
            overlap_start = max(var['start'], ref_start)
            overlap_end   = min(var['end'], ref_end)
            if overlap_start < overlap_end:
                overlap_len = overlap_end - overlap_start
                de_len = var['end'] - var['start']
                ref_len = ref_end - ref_start
                overlap_pct = (overlap_len / float(max(de_len, ref_len))) * 100.0
                if overlap_pct >= overlap_threshold:
                    is_found = True
                    break
        if is_found:
            found.append(var)
        else:
            not_found.append(var)
    return found, not_found

def init_worker(ref_data):
    global ref_trees
    ref_trees = ref_data

def process_file(file_path):
    de_novo_vars = read_variants(file_path)
    found, not_found = find_overlaps_in_reference(de_novo_vars, ref_trees)
    return (file_path, found, not_found)

# Directory with de novo VCF files and reference VCF file
denovo_dir = "/data/denovo/both_delly_lumpy_denovo/"
gnomad_vcf_file = "/data/denovo/gnomad.v4.1.sv.sites.vcf"
output_dir = "/data/denovo/new/"
found_vcf_path = os.path.join(output_dir, "de_novo_found_in_gnomad.vcf")
not_found_vcf_path = os.path.join(output_dir, "de_novo_not_found_in_gnomad.vcf")
os.makedirs(output_dir, exist_ok=True)

# Clean previous output files
for f in [found_vcf_path, not_found_vcf_path]:
    if os.path.exists(f):
        os.remove(f)

print("Loading gnomAD variants...")
gnomad_variants = read_variants(gnomad_vcf_file)
print(f"Loaded {len(gnomad_variants)} gnomAD variants.")
print("Building interval trees...")
reference_trees = build_interval_trees_by_chrom(gnomad_variants)
print(f"Built interval trees for {len(reference_trees)} chromosomes.")

file_list = glob.glob(os.path.join(denovo_dir, "*.vcf"))
if not file_list:
    print("No VCF files found.")
    exit()

cpu_count = min(multiprocessing.cpu_count(), len(file_list), 8)
print(f"Starting processing on {cpu_count} parallel workers...")
with multiprocessing.Pool(processes=cpu_count, initializer=init_worker, initargs=(reference_trees,)) as pool:
    results = pool.map(process_file, file_list)

# Collect variant records
total_found = 0
total_not_found = 0
total_vcfs = len(results)
found_records_by_header = []
not_found_records_by_header = []

for file_path, found_list, not_found_list in results:
    found_keys = {(v['chrom'], v['start'], v['end']) for v in found_list}
    not_found_keys = {(v['chrom'], v['start'], v['end']) for v in not_found_list}

    with pysam.VariantFile(file_path) as in_vcf:
        header = in_vcf.header.copy()
        for rec in in_vcf:
            key = (rec.chrom, rec.start, rec.stop)
            if key in found_keys:
                found_records_by_header.append((header, rec))
                total_found += 1
            elif key in not_found_keys:
                not_found_records_by_header.append((header, rec))
                total_not_found += 1

# Write output
if found_records_by_header:
    with pysam.VariantFile(found_vcf_path, 'w', header=found_records_by_header[0][0]) as found_out:
        for header, rec in found_records_by_header:
            try:
                found_out.write(rec)
            except Exception as e:
                print(f"Skipping record due to write error: {e}")

if not_found_records_by_header:
    with pysam.VariantFile(not_found_vcf_path, 'w', header=not_found_records_by_header[0][0]) as not_found_out:
        for header, rec in not_found_records_by_header:
            try:
                not_found_out.write(rec)
            except Exception as e:
                print(f"Skipping record due to write error: {e}")

# Report
print("\n=== Final Report ===")
print(f"Total VCF files processed: {total_vcfs}")
total_variants = total_found + total_not_found
if total_variants > 0:
    print(f"Total variants processed: {total_variants}")
    print(f"Variants found in gnomAD: {total_found} ({(total_found / total_variants) * 100:.2f}%)")
    print(f"Variants not found in gnomAD: {total_not_found} ({(total_not_found / total_variants) * 100:.2f}%)")
else:
    print("No variants were processed. Check if the gnomAD file contains valid structural variants.")