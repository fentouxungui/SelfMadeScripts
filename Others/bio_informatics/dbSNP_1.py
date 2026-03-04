import xlrd
import xlwt
import vcf
import time

vcf_reader = vcf.Reader(open('C:\\Users\Zhang Yongchao\Desktop\Test\dbsnp_v150_hg19_refseq.vcf', 'r'))
workbook=xlrd.open_workbook('C:\\Users\Zhang Yongchao\Desktop\Test\SNP+call+file.xlsx') ## 打开文件
sheet=workbook.sheet_by_index(0)# 根据sheet索引或者名称获取sheet内容
wb_new = xlwt.Workbook()
ws_new = wb_new.add_sheet('A Test Sheet',cell_overwrite_ok=True)
dictA={}# 创建一个字典，用于查找
dictB={}# 创建一个字典，用于查找
for i in range(1,len(sheet.col_values(1))):
    dictA[sheet.row_values(i)[1]]=i #关联字典的键和值

for record in vcf_reader:
    bb = record.ID.split(';')
    for b in bb:
        dictB[b]=record.POS

for j in dictA:
    if j in dictB: #如果当前读取的rs号在已有的字典中存在，则输出
        ws_new.write(j, 0, record.POS)

wb_new.save('example.xls')

print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
