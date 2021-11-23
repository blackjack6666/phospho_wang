from selenium import webdriver
from selenium.webdriver.support.ui import Select
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
import time

driver = webdriver.Chrome(executable_path='C:/tools/chromedriver.exe')
url = 'https://uic.ca1.qualtrics.com/jfe/form/SV_2fem8PpxGd7Ab54'
driver.get(url)
driver.maximize_window()

wait = WebDriverWait(driver, 5)
x_path = '//*[@id="QR~QID1~1"]/option[12]'

select = Select(driver.find_element(By.ID,'QR~QID1~1'))
select.select_by_value('17')

time.sleep(1)

select2 = Select(driver.find_element(By.ID,'QR~QID2~1'))
select2.select_by_value('11')

time.sleep(1.5)

select3 = Select(driver.find_element(By.ID,'QR~QID3~1'))
select3.select_by_value('1')

time.sleep(3)
### click by id
driver.find_element(By.ID,'NextButton').click()

