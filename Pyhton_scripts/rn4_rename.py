## File to add gene names to GTF files from table browser
import fileinput
import sys
import re



#opening file with transcripts IDs to gene names
annotation = open("/home/glk10/Desktop/rn4_refseq_Name-ID.txt").readlines()
#create dictionary of transcript ID: gene names
gtf = {}
for line in annotation:
    x = line.split("\t")
    x[1] = x[1].strip()
    gtf[x[0]] = [x[1]]

#replacing transcript ID (in gene name position) with gene names
with open("rn4_refseq_rename.gtf", "w") as new_file:
    with open("rn4.gtf") as old_file:
        count = 0
        for line in old_file:
                m = re.search(r'gene_id "(\w+)";',line)
                id = m.group(1)
                if m:
                    replace = str(gtf[m.group(1)])
                    split = replace.split("'")
                    og = 'gene_id "{}"'.format(id)
                    new = 'gene_id "{}"'.format(split[1])
                    new_file.write(line.replace(og, new))
                    count = count + 1
                    print(count)
        sys.stdout.write(line)
