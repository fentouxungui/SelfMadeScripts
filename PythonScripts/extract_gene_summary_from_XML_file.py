# ref: https://blog.51cto.com/u_14782715/5082276
from Bio import Entrez
# =====解析大文件=====
# hd_parse = open("C:\\Users\Xi_Lab\\Desktop\\NCBI-extract-Summary\\Homo_sapiens.xml", "rb")
hd_parse = open("C:\\Users\Xi_Lab\\Desktop\\NCBI-extract-Summary\\subset.xml", "rb")
res_parse = Entrez.parse(hd_parse)
for record in res_parse:
    status = record['Entrezgene_track-info']['Gene-track']['Gene-track_status']
    if status.attributes['value'] == 'discontinued':
        continue
    gene_id = record['Entrezgene_track-info']['Gene-track']['Gene-track_geneid']
    # 有一些ID没有Gene-ref_locus
    if 'Gene-ref_locus' not in record['Entrezgene_gene']['Gene-ref'].keys():
        continue
    else:
        gene_name = record['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']
    # 有的没有Entrezgene_summary
    if 'Entrezgene_summary' in record.keys():
        gene_summary = record['Entrezgene_summary']
    else:
        gene_summary = ""
    with open('C:\\Users\Xi_Lab\Desktop\\NCBI-extract-Summary\\NCBI_gene_summary_2.txt', 'a') as f:
        f.write(gene_id + '\t' + gene_name + '\t' + gene_summary + '\n')
    # print(gene_id, gene_name, gene_summary)

# 跑到一个地方报错了，可以提取剩余的
# 获取截取位置
# cat Homo_sapiens.xml | grep "<Entrezgene>\|<Gene-ref_locus>" -n | grep "LRRC56" -B10 -A10
# cat Homo_sapiens.xml | sed '4,425058435d' > subset.xml

# 想到另外一种提取办法
# cat Homo_sapiens.xml | grep "<Gene-track_geneid>\|<Gene-ref_locus>\|<Entrezgene_summary>" > subset.txt

