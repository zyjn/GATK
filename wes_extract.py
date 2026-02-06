import os
import subprocess

input_dir = "/mnt/project/Bulk/DRAGEN WGS/DRAGEN population level WGS variants, PLINK format [500k release]/" 
plink_exe = "/opt/notebooks/plink2"
gene_list_file = "gene_list.txt"
output_dir = "./Ensembl_VCFs_DRAGEN/"
file_prefix_pattern = "ukb24308_c{chr}_b0_v1" # DRAGEN的前缀通常是 ukb24308
master_psam = "./fixed_samples.psam" 

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if not os.path.exists(master_psam):
    print(f"{master_psam}")
    exit(1)

# read.list
genes_by_chr = {}
with open(gene_list_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 5:
            c = parts[0].replace('chr', '').replace('X','X').replace('Y','Y')
            if c not in genes_by_chr: genes_by_chr[c] = []
            genes_by_chr[c].append({
                'start': parts[1], 'end': parts[2], 'id': parts[4], 'symbol': parts[3]
            })

sorted_chroms = sorted(genes_by_chr.keys(), key=lambda x: int(x) if x.isdigit() else 999)

for chrom in sorted_chroms:
    print(f"processing: chromosome {chrom}...") 
    pfile_base = os.path.join(input_dir, file_prefix_pattern.format(chr=chrom))
    pgen_file = pfile_base + ".pgen"
    pvar_file = pfile_base + ".pvar"
    
    if not os.path.exists(pgen_file):
        print(f"can't find: {pgen_file}")
        continue
    if not os.path.exists(pvar_file):
        print(f"can't find: {pvar_file}")
        continue

    count = 0
    for gene in genes_by_chr[chrom]:
        out_name = os.path.join(output_dir, gene['id'])
        
        if os.path.exists(out_name + ".vcf.gz") and os.path.getsize(out_name + ".vcf.gz") > 100:
            count += 1; continue       
        cmd = [
            plink_exe,
            "--pgen", pgen_file,       # .pgen
            "--pvar", pvar_file,       # .pvar
            "--psam", master_psam,     # .psam
            "--chr", chrom,
            "--from-bp", gene['start'],
            "--to-bp", gene['end'],
            "--export", "vcf", "bgz",
            "--out", out_name,
            "--silent",
            "--memory", "12000"
        ]
        
        try:
            # print("Running:", " ".join(cmd)) 
            subprocess.run(cmd, check=True)
            count += 1
        except subprocess.CalledProcessError as e:
            print(f"---fail: {gene['symbol']} (ID: {gene['id']})")
            # 如果出错，打印一下出错的文件路径，确认没找错文件
            print(f"test PGEN: {pgen_file}")

    print(f"chromosome {chrom} success: {count}/{len(genes_by_chr[chrom])}")

print("done")