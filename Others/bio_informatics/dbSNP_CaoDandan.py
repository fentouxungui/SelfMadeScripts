import xlrd
import xlwt

workbook=xlrd.open_workbook('C:\\Users\Zhang Yongchao\Desktop\Test\SNP+call+file.xlsx') ## 打开文件
print(workbook.sheet_names()) # [u'sheet1', u'sheet2']   获取所有sheet
sheet=workbook.sheet_by_index(0)# 根据sheet索引或者名称获取sheet内容
print(sheet.name,sheet.nrows,sheet.ncols) # sheet的名称，行数，列数

# 获取整行和整列的值（数组）
rows = sheet.row_values(3) # 获取第四行内容
cols = sheet.col_values(2) # 获取第三列内容
print(rows)
print(cols)
# 获取单元格内容
print(sheet.cell(1, 0).value.encode('utf-8'))
print(sheet.cell_value(1, 0).encode('utf-8'))
print(sheet.row(1)[0].value.encode('utf-8'))
# 获取单元格内容的数据类型   ctype : 0 empty,1 string, 2 number, 3 date, 4 boolean, 5 error
print(sheet.cell(1,0).ctype)


dict={}# 创建一个字典，用于查找
for i in range(1,len(sheet.col_values(1))):
    dict[sheet.col_values(1)[i]]=sheet.row_values(i) #关联字典的键和值
f=open('C:\\Users\Zhang Yongchao\Desktop\input.txt','r')#打开索引的文件
b=open('C:\\Users\Zhang Yongchao\Desktop\output.txt','a')#打开输出文件
lines = f.readlines() #读取dbsnp中的所有行
for line in lines:  #逐行读取
    a=line.split()  #分解每一行获得独立元素
    if a[0] in dict:    #如果当前读取的rs号在已有的字典中存在，则输出
        print(dict[a[0]][0],dict[a[0]][1],dict[a[0]][2],a[1],a[2],file=b)
f.close()
b.close()