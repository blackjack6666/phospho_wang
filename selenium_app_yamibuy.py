from selenium import webdriver
from selenium.webdriver.support.ui import Select
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.support import expected_conditions as EC
import time
import requests
import tweepy
from twitter_keys import *
from twython import Twython


def lauch_browser(url):
    ser = Service('C:/tools/chromedriver.exe')
    op = webdriver.ChromeOptions()
    driver = webdriver.Chrome(service=ser, options=op)
    driver.get(url)
    driver.maximize_window()
    return driver


def login(email,password):
    # email input x path
    email_ele = driver.find_element(By.XPATH,'//*[@id="login-container"]/div[1]/div[2]/input')
    email_ele.clear()
    email_ele.send_keys(email)
    time.sleep(1)
    # password input x path
    password_ele = driver.find_element(By.XPATH,'//*[@id="login-container"]/div[1]/div[3]/input')
    password_ele.clear()
    password_ele.send_keys(password)


def send_twitter(message):
    client = tweepy.Client(bearer_token=bearer_token, consumer_key=consumer_key, consumer_secret=consumer_secret,
                           access_token=access_token, access_token_secret=access_token_secret)
    # api = tweepy.API(client)
    client.create_tweet(text=message)

    # recipient_id = 'beefnoodles111'
    # api.send_direct_message(recipient_id, message)
    print (f'sending tweets: {message}')


if __name__=='__main__':

    # yamibuy item URL

    url = 'https://www.yamibuy.com/zh/p/mini-la-tiao-soybean-snacks-360g/1017041441?track=rec-branchpage-58-1017041441&scene=cms_channel.snack&bu_type=display&module_name=cms_component_items&index=58&content=1017041441'
    driver = lauch_browser(url)
    item = url.split('content=')[1]

    quehuo_Xpath = "//button[@class='action-add_cart float_left button-solid-yellow-active']"
    que_huo = True

    # keep checking if item is missing or not, until it's in stock
    while True:
        try:
            driver.find_element(By.XPATH,quehuo_Xpath)

            print ('que huo...')
            time.sleep(86400) # wait a day and check again
        except:
            print ('bu que huo, go on auto-buying...')
            que_huo = False
            break

    if not que_huo:  # if itme is in stock
        x_path = '//*[@id="item-container"]/div[3]/div[2]/div[4]/div[1]/div[4]/div[1]/i[2]' # add item number

        for i in range(5): # buy 5 times
            time.sleep(1)

            WebDriverWait(driver, 5).until(EC.element_to_be_clickable((By.XPATH, x_path))).click() # add 1 every 3 seconds



        #add to cart
        add_to_cart_x_path = '//*[@id="item-container"]/div[3]/div[2]/div[4]/div/div[4]/div[2]/button'
        time.sleep(1)
        WebDriverWait(driver, 10).until(EC.element_to_be_clickable((By.XPATH, add_to_cart_x_path))).click()

        #close poped up window
        window_close_xpath = '//*[@id="item-container"]/div[7]/span/i'
        time.sleep(2)
        WebDriverWait(driver, 3).until(EC.element_to_be_clickable((By.XPATH, window_close_xpath))).click()


        # go to checkout window
        check_out_xpath = '//*[@id="header-cart"]/a/span'
        time.sleep(2)
        WebDriverWait(driver, 3).until(EC.element_to_be_clickable((By.XPATH, check_out_xpath))).click()

        # checkout
        check_out_xpath2 = '//*[@id="cart-container"]/div[6]/a'
        time.sleep(2)
        WebDriverWait(driver, 3).until(EC.element_to_be_clickable((By.XPATH, check_out_xpath2))).click()

        # login with username and password
        email,password = 'sxhlzcn@qq.com','sxh8836995...'
        login(email,password)

        # send twitter
        send_twitter('%s added cart and logged in' % item)  # it will show 403 forbidden if tweet sent is same


    driver.close()