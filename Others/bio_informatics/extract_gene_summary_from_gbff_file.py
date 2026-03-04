import os
def Merge(dict1, dict2):
    res = {**dict1, **dict2}
    return res
def extract_summary(file):
    f = open(file)
    locus2comment = {}
    in_comment = False
    for line in f:
        if line[0:5] == "LOCUS":
            locus = line.split()[1]
            comment = ""
        elif line[0:7] == "COMMENT":
            in_comment = True
            comment += line.split("    ")[1].replace("\n", " ")
        elif line[0:7] == "PRIMARY":
            in_comment = False
            try:
                locus2comment[locus] = comment.split("Summary:")[1]
            except:
                locus2comment[locus] = comment
        elif in_comment:
            comment += line.split("            ")[1].replace("\n", " ")
    return locus2comment

dirpath = 'C:\\Users\Xi_Lab\Desktop\\NCBI-extract-Summary\DownloadedData'

res = {}
for f in os.listdir(dirpath):
    if os.path.splitext(f)[-1] == ".gbff":
        print("Processing: " + f + "...")
        res = Merge(res, extract_summary(file=dirpath + "\\" + f))

for locus in sorted(res):
    with open('C:\\Users\Xi_Lab\Desktop\\NCBI-extract-Summary\\NCBI_gene_summary.txt', 'a') as f:
        f.write(locus + '\t' + res[locus] + '\n')
    # print(locus + '\t' + res[locus])