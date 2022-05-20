from urllib.request import Request, urlopen
from selenium_app_yamibuy import send_twitter
import time
import random
from datetime import datetime


while True:
    req = Request('https://www.jellycat.com/us/munro-scottie-dog-mun3sd/',
                  headers={'User-Agent': 'Mozilla/5.0'})
    webpage = urlopen(req).read().decode('utf-8')
    item = webpage.find('Out Of Stock')
    if item != -1:
        print(f'sleeping...current time: {datetime.now()}')
        time.sleep(3600+random.randint(1,20))

    else:
        send_twitter('scottie dog jellycat available https://www.jellycat.com/us/munro-scottie-dog-mun3sd/')
        break
print ('loop break, twitter sent')

