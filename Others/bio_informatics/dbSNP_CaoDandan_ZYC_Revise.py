#CPU i5 四核 20%左右，RAM300-500M左右
import xlrd
import xlwt
import vcf
import time

print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))

vcf_reader = vcf.Reader(open('C:\\Users\Zhang Yongchao\Desktop\Test\dbsnp_v150_hg19_refseq.vcf', 'r'))
workbook=xlrd.open_workbook('C:\\Users\Zhang Yongchao\Desktop\Test\SNP+call+file.xlsx') ## 打开文件
sheet=workbook.sheet_by_index(0)# 根据sheet索引或者名称获取sheet内容
wb_new = xlwt.Workbook()
ws_new = wb_new.add_sheet('A Test Sheet',cell_overwrite_ok=True)
dict={}# 创建一个字典，用于查找
for i in range(1,len(sheet.col_values(1))):
    dict[sheet.row_values(i)[1].strip()]=i #关联字典的键rsID(去除空格)和值（行坐标）

for record in vcf_reader:
    a=record.ID.split(';')
    for b in a:
        if b in dict: #如果当前读取的rs号,并去除空格，在已有的字典中存在，则输出
            ws_new.write(dict[b], 0, record.CHROM)
            ws_new.write(dict[b], 1, record.POS)

wb_new.save('example.xls')

print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))