#!/usr/bin/env python
# coding:utf-8
import re
from collections import Counter
FILESOURCE = "C:\\Users\\Zhang Yongchao\\Desktop\\insert.mm.dist.among.unique.scaled.mOTU"
def getMostCommonWord(articlefilesource):
	'''输入一个英文的传文本文件，统计其中的单词出现的个数'''
	pattern = r'''COG00\d'''
	with open(articlefilesource) as f:
		r = re.findall(pattern,f.read())
		return Counter(r)
if __name__ == '__main__':
	print(getMostCommonWord(FILESOURCE))