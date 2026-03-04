import xlrd

workbook = xlrd.open_workbook('C:\\Users\Zhang Yongchao\Desktop\Test\SNP+call+file.xlsx')
sheet=workbook.sheet_by_index(0)
dict={}# 创建一个字典，用于查找

for i in range(1,len(sheet.col_values(1))):
    dict[sheet.col_values(1)[i].strip()]=sheet.row_values(i) #关联字典的键和值


b=open('C:\\Users\Zhang Yongchao\Desktop\SNP call file copy.txt','a')#打开输出文件
d=open('C:\\Users\Zhang Yongchao\Desktop\Test\dbsnp_v150_hg19_refseq.txt','r')#打开索引的文件

lines = d.readlines() #读取dbsnp中的所有行
for line in lines:  #逐行读取
    if '#' not in line:
        a=line.split()
        aa=a[2].split(';')
        for aaa in aa:#分解每一行获得独立元iiii
            if aaa in dict:    #如果当前读取的rs号在已有的字典中存在，则输出
                print(dict[aaa][0],dict[aaa][1],dict[aaa][2],a[0],a[1],file=b)
d.close()
b.close()
