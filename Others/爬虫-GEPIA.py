from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.select import Select
import pandas as pd
from selenium.webdriver.support.ui import WebDriverWait

# driver = webdriver.Firefox()
driver = webdriver.Chrome(executable_path="C:\\Program Files\\Git\mingw64\\bin\\chromedriver.exe")
driver.get("http://gepia.cancer-pku.cn/detail.php")
driver.implicitly_wait(10) # implicitly_wait()隐式等待网页加载完成

driver.find_element_by_id("a_degenes").click() # 单击差异分析选项卡

driver.find_element_by_id("degenes_limma").click() # 设置差异分析方法为limma

select = driver.find_element_by_id("degenes_overlow")
Select(select).select_by_index(0) # 选择染色体分布为 过表达

select = driver.find_element_by_id("degenes_select")
Select(select).select_by_index(23) # 选择数据集 5:COAD; 23:READ

driver.find_element_by_name("degenes_submit").click() # 提交

select = driver.find_element_by_name("table_id_example_length")
Select(select).select_by_index(3) # 设置每页显示100个


webtable_df = pd.read_html(driver.find_element_by_id("table_id_example").get_attribute('outerHTML'))[0]
# webtable_df.to_csv('file2.csv')

while (driver.find_element_by_link_text('Next').get_attribute("class") == 'paginate_button next'):
    driver.find_element_by_id('table_id_example_next').click()
    webtable_df = pd.concat([webtable_df, pd.read_html(driver.find_element_by_id("table_id_example").get_attribute('outerHTML'))[0]], axis = 0)
webtable_df.to_csv('READ.up.genes.csv')

driver.close()