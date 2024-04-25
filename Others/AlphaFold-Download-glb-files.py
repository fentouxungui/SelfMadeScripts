from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.select import Select
import pandas as pd
import time
import os
import time
# os.chdir("C:\\Users\\Xi_Lab\Desktop\\neuropeptide\\果蝇GutEEs\\Xi_Filtered")
os.chdir("C:\\Users\\Xi_Lab\\Desktop\\neuropeptide\\人GutEECs\\数据整合")
# driver = webdriver.Firefox()
# driver = webdriver.Chrome(executable_path="C:\\Program Files\\Git\mingw64\\bin\\chromedriver.exe")
s = Service(r"C:\\Program Files\\Git\mingw64\\bin\\chromedriver.exe")
driver = webdriver.Chrome(service=s)

protein_df = pd.read_csv("Neuropeptide-Proteins.csv")
proteins = protein_df["Protein"]
def get_files_with_postfix(folder_path, postfix):
    files = []
    for file_name in os.listdir(folder_path):
        if file_name.endswith(postfix):
            files.append(file_name)
    return files


for i in proteins:
    if i in get_files_with_postfix("C:\\Users\\Xi_Lab\\Downloads",postfix=".glb"):
        continue
    print(i)
    driver.get("https://alphafold.ebi.ac.uk/entry/" + i)
    driver.implicitly_wait(10)  # implicitly_wait()隐式等待网页加载完成
    # driver.find_element(By.XPATH, '//*[@title="glTF 2.0 Binary (.glb)"]').click()
    # driver.find_element(By.XPATH, '//button[@title="glTF 2.0 Binary (.glb)"]').click()
    # driver.find_element(By.XPATH, '//div[@class="msp-layout-region msp-layout-right"]//button[@class="msp-btn msp-btn-block msp-no-overflow msp-action-menu-button"]').click()
    ExportCollapsed = True
    # 有的时候“driver.find_element(By.XPATH, '//button[text()="Export Geometry"]').click()”代码运行不成功，需要重复运行几次，或者用鼠标操作，定位到那个区域，不用点就可以！奇怪~
    while ExportCollapsed:
        try:
            driver.find_element(By.XPATH, '//button[text()="Export Geometry"]').click()
            ExportCollapsed = False
        except:
            time.sleep(2)
    driver.find_element(By.XPATH, '//button[text()="Save"]').click()
time.sleep(3)
driver.quit()
