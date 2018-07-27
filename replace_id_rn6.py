import fileinput
import sys
import re
from itertools import islice

# Create a dictionary wich include paired StringTie id and RefSeq id.
refmap = open("merged.stringtie_merged_rn6.gtf.refmap")

ids = {}
for line in islice(refmap, 1, None):
    
    d = line.split("\t")
    d[3] = d[3].strip().split("|")
    stid = str(d[3][:-1]).replace("'","").replace("[","").replace("]","")
    gname = str(d[0])
    ids[stid]= gname


# Write a new result file replace StringTie id into RefSeq id.
with open("rn6_geneid_result.csv", "w") as new_file:
    with open("rn6_genes_results.csv") as old_file:
        new_file.write('"feature","id","fc","pval","qval"'+'\n')
        count=0
        for line in islice(old_file,1,None):
            splitline = line.split(",")
            oldid = splitline[1].replace('"',"")
            if oldid in ids:
                new= ids.get(oldid)
                new_file.write(line.replace(oldid,new))
                count = count+1
                print(count)

               
                
            


                        

