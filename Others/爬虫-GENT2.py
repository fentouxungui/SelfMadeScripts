# GENT2: An Updated Gene Expression Database for Normal and Tumor Tissues

# Introduction

# Gene Expression database of Normal and Tumor tissues 2 (GENT2) is an updated version of GENT, which has provided a
# # user-friendly search platform for gene expression patterns across different normal and tumor tissues compiled from
# public gene expression data sets. We refactored GENT2 with recent technologies such as Apache Lucene indexing for fast
# search and Google Web Toolkit (GWT) framework for a user-friendly web interface. Now, GENT2 contains more than 68,000
# samples and has several new useful functions. First, GENT2 now provides gene expression across 72 different tissues
# compared to 57 in GENT. Second, with increasing importance of tumor subtypes, GENT2 provides an option to study the
# differential expression and its prognostic significance based on tumor subtypes. Third, whenever available, GENT2
# provides prognostic information of a gene of interest. Fourth, GENT2 provides a meta-analysis of survival information
# to provide users more reliable prognostic value of a gene of interest. In conclusion, with these significant
# improvements, GENT2 will continue to be a useful tool to a wide range of researchers. GENT2 is freely available at
# http://gent2.appex.kr.

# Keywords: Cancer subtype profiling; Large scale microarray web-database; Meta-survival analysis; Survival analysis;
# Tissue and cell line wide gene expression profiling.

# Please cite your use of GENT2 in your publication:

# Park SJ, Yoon BH, Kim SK*, Kim SY*. GENT2: an updated gene expression database for normal and tumor tissues.
# BMC Med Genomics. 2019 Jul 11;12(Suppl 5):101. doi: 10.1186/s12920-019-0514-7.

# 注意： GENE2 提供了全部数据的下载链接，但是文件太大，约20GB，我们仅需检索几十个TFs的信息，所以在这里用爬虫方式进行提取。

from selenium import webdriver
from selenium.webdriver.common.by import By


driver = webdriver.Chrome(executable_path="C:\\Program Files\\Git\mingw64\\bin\\chromedriver.exe")
driver.get("http://gent2.appex.kr/gent2/")
driver.implicitly_wait(10) # implicitly_wait()隐式等待网页加载完成

# 点击 'Search -Gene profile'
driver.find_element_by_id("x-auto-13-input").send_keys("TK")

driver.find_element_by_class_name("GD3CNKB-j-c").click()
driver.find_element_by_class_name("GD3CNKB-f-d x-toolbar x-small-editor").click()
driver.find_element_by_class_name("GD3CNKB-t-b").click()
driver.find_element_by_class_name("GD3CNKB-u-j GD3CNKB-u-d GD3CNKB-u-q GD3CNKB-u-n").click()
driver.find_element_by_class_name("GD3CNKB-r-c").click()
driver.find_element_by_class_name("GD3CNKB-r-d GD3CNKB-r-g  GD3CNKB-r-o").click()
driver.find_element_by_class_name("GD3CNKB-r-g").click()

driver.execute_script('this.__gwtLastUnhandledEvent="load"')
#driver.execute_script("alert('hello')")