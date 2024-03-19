# 此代码未曾测试
# https://gist.github.com/mmendez12/8b87f0112f6e203d4c81
from Bio import Entrez
def human_gene_summary_from_entrez(gene_name='RRP9'):
    Entrez.email = "toto@tata.com" # Always tell NCBI who you are
    term = "{}[Gene Name] AND Homo Sapiens[Organism]".format(gene_name)
    handle = Entrez.esearch(db='gene', term=term)
    result = Entrez.read(handle)
    handle.close()

    gene_id = result["IdList"][0]

    handle = Entrez.esummary(db="gene", id=gene_id)
    result = Entrez.read(handle)
    handle.close()

    res = result['DocumentSummarySet']['DocumentSummary'][0]
    return res['Name'], res['Summary']