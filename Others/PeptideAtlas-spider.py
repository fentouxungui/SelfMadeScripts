from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.select import Select
import pandas as pd
import time
import os

os.chdir("C:\\Users\\Xi_Lab\Desktop\\neuropeptide\\人GutEECs\\scRNAseq\\PublicData-27")
# driver = webdriver.Firefox()
# driver = webdriver.Chrome(executable_path="C:\\Program Files\\Git\mingw64\\bin\\chromedriver.exe")
s = Service(r"C:\\Program Files\\Git\mingw64\\bin\\chromedriver.exe")
driver = webdriver.Chrome(service=s)
driver.get("https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProtein")
db_species = Select(driver.find_element(By.NAME, 'atlas_build_id'))  # 实例化Select
db_species.select_by_value("550")

protein_df = pd.read_csv("highly-expressed-genes-2-Annotated.csv")
proteins_raw = protein_df["HGNC_UniProt.ID.supplied.by.UniProt"].dropna().values
proteins = []

peptideatlas_db_local = "C:\\Users\\Xi_Lab\Desktop\\neuropeptide\\人GutEECs\\PeptideAtlas\\Spider-Observed-in-Experiments"
proteins_local = [i.replace(".csv", "") for i in os.listdir(peptideatlas_db_local)]

for a_proteins_raw in proteins_raw:
    proteins = proteins + a_proteins_raw.split(", ")

print("Total found ", len(proteins), " proteins.")
print(len(list(set(proteins).intersection(set(proteins_local)))), " proteins are found in local.")
proteins = list(set(proteins) - set(proteins_local))
print(len(proteins), " proteins will be downloaded.")

start = time.perf_counter()
scale = len(proteins)
proteins_not_found = []

for i in range(scale + 1):
    # 一个*或-表示一个样本！
    a = "*" * i
    b = "." * (scale - i)
    c = (i / scale) * 100
    dur = time.perf_counter() - start
    print("\r{:^3.0f}%[{}->{}]{:.2f}s".format(c, a, b, dur), end="")
    # ------------ function start
    protein = proteins[i - 1]
    driver.find_element(By.NAME, "protein_name").clear()
    driver.find_element(By.NAME, "protein_name").send_keys(protein)
    driver.find_element(By.NAME, "action").click()
    driver.implicitly_wait(10)  # implicitly_wait()隐式等待网页加载完成
    try:
        driver.find_element(By.ID, "getprotein_overview_div")
        try:
            driver.find_element(By.LINK_TEXT, '[Show more rows]').click()
        except:
            pass
        try:
            webtable_df = pd.read_html(
                driver.find_element(By.XPATH, "//div[@id='getprotein_samplelist_div']//table").get_attribute(
                    'outerHTML'))[0]
            webtable_df = webtable_df.drop(webtable_df.index[[0, 1, 2]])
            # webtable_df.columns.values
            webtable_df = webtable_df.rename(columns={'Experiment ID\xa0▿': 'Experiment ID'})
            webtable_df.to_csv(peptideatlas_db_local + '\\' + protein + '.csv')
        except:
            proteins_not_found.append(protein)
    except:
        proteins_not_found.append(protein)

driver.quit()
if len(proteins_not_found) != 0:
    print("\nTotal ", len(proteins_not_found), " proteins are not found. they are:")
    print(proteins_not_found)
