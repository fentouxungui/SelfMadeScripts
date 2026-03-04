import pip
import xlwt
import xlrd
import vcf #导入PyVCF库，这个数据库是下载的一系列Gene信息vcf格式文件
#调用Reader对象处理vcf文件
vcf_reader = vcf.Reader(open('C:\\Users\Zhang Yongchao\Desktop\Test\dbsnp_v150_hg19_refseq.vcf', 'r'))
workbook = xlrd.open_workbook('C:\\Users\Zhang Yongchao\Desktop\Test\SNP+call+file.xlsx', 'rb')
table = workbook.sheet_by_name('Sheet1')
print(vcf_reader.fetch(1))


#新建一个Excel
wb_new = xlwt.Workbook()
ws_new = wb_new.add_sheet('A Test Sheet')
# 逐条获取查询结果中的数据，返回的是一个Record(Gene=a, ID=b, POS=c)

records = next(vcf_reader)
for record in records:    #print(record.ID)
    #print(type(record.ID))
    for j in range(1, table.nrows-2):
        #print(table.cell(j,1).value)
        if record.ID == table.cell(j,1).value:
        # 第i行第3列写入内容
            ws_new.write(j, 0, record.POS)
            print(record.POS)
            break
wb_new.save('example.xls')
