import xlrd
import xlwt
import vcf
import time

vcf_reader = vcf.Reader(open('C:\\Users\Zhang Yongchao\Desktop\Test\dbsnp_v150_hg19_refseq.vcf', 'r'))
print(vcf_reader)
print(next(vcf_reader))