#!/bin/python3

import re
from urllib.request import urlopen
import os


# find the FTP address from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GEO
response = urlopen("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81916")
pattern = re.compile("<a href=\"(.*?)\">\(ftp\)</a>")
# use wget from shell to download SRA data
ftp_address = re.search(pattern, response.read().decode('utf-8'))
print(ftp_address)
os.system(' wget -nd -r 1 -A *.sra ' + ftp_address)
'''
"""
用sys.argv从命令行中读取参数
用urllib.request向网页发起请求，获取response
用正则表达式(re)提取FTP地址
用os.system运行shell的命令
"""

