#!/usr/bin/env python



#from hidden.py import password
from selenium import webdriver
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.support import expected_conditions as EC
import time


#driver = webdriver.Chrome('/Users/jeremyjacobson/.wdm/drivers/chromedriver/mac64/96.0.4664.45/chromedriver')
driverService = Service('/Users/jeremyjacobson/.wdm/drivers/chromedriver/mac64/96.0.4664.45/chromedriver')   
driver = webdriver.Chrome(service=driverService)         


#driver.get('https://login.alterramtnco.com/login?state=hKFo2SB1Q0p1Y25IdW1pVkZMVTVVU2h2QXpySXl5YkRack5vSqFupWxvZ2luo3RpZNkgSTRJVGJQeE9kd3lWWDhURTJoY2t0S29vd2JaUWRJQjKjY2lk2SBvQ2kwRlFiRDR5SVhldWdGWU43T3NBd2cwckFINDAzVg&client=oCi0FQbD4yIXeugFYN7OsAwg0rAH403V&protocol=oauth2&code_challenge=fGceKVaFLdwV2j6PTgC98aY_55vU9quf3IsNNrG0aQI&code_challenge_method=S256&leeway=60&nonce=458e1dc5f46669e2ed9df5c51600f14d&redirect_uri=https%3A%2F%2Faccount.ikonpass.com%2Fauth%2Fauth0%2Fcallback%3Fui_locales%3Den%26redirect_uri%3Dhttps%3A%2F%2Fikonpass.com%2F&response_type=code&scope=openid%20profile&ui_locales=en#/')
driver.get("https://account.ikonpass.com/en/myaccount/add-reservations/")
driver.implicitly_wait(15)

#Click the Log In button
#driver.find_element(By.CLASS_NAME, "inner").click

#saved locations
Jackson = "sc-prRrb eYneqP"
Taos = "sc-prRrb ibDMsT"

test_date = "Thu Jan 20 2022"



myPageTitle = driver.title              
print(myPageTitle) 
#assert "Ikon" in myPageTitle   

my_email = "Jeremy.Jacobson3402@gmail.com"
my_password = password

email_location = driver.find_element(By.NAME, "email")

password_location = driver.find_element(By.NAME, "password")

email_location.send_keys(my_email)
password_location.send_keys(my_password)


# Log in
driver.find_element(By.CLASS_NAME, "inner").click()

#Navigate to Registration Page - opens on Jan 4th.
#driver.get("https://www.ikonpass.com/en/destinations/crystal-mountain")

#driver.get("https://account.ikonpass.com/en/myaccount/add-reservations/")
time.sleep(10)


#select location
driver.find_element(By.CLASS_NAME, Jackson).click()
driver.find_element(By.CLASS_NAME, "sc-AxjAm jxPclZ sc-qQyIJ bHttEQ").click()

#date
date = "\"[aria-label=" + test_date + "]\""
driver.find_elements_by_css_selector(date)


#If Dismiss button is present, repeat in 91 seconds.
Res_availibility = False
Res_availibility = driver.find_element(By.CLASS_NAME, "sc-AxjAm jxPclZ sc-pYdvZ uunor")
attempts = 1
while Res_availibility == False:
    time.sleep(5)
    print("Attempts to reserve: ", attempts)


#save buttom
driver.find_element(By.CLASS_NAME, "sc-AxjAm jxPclZ sc-pleUP jkQaba").click()

#First Confirmation button
driver.find_element(By.CLASS_NAME, "sc-AxjAm jxPclZ sc-qOhCc ijLuHQ").click()

#Checkbox
driver.find_element(By.CLASS_NAME, "sc-qWdWY dCFeja").click()

#Final Confirmation button
driver.find_element(By.CLASS_NAME, "sc-AxjAm jxPclZ").click()





#Final web page - could use for confirmation
#https://account.ikonpass.com/en/myaccount/reservation-success



#Menu Drop Down   
# driver.find_element(By.CLASS_NAME, "navigation-menu-icon").click()
# driver.find_element(By.CLASS_NAME, "active amp-button unstyled").click()
# driver.find_element(By.CLASS_NAME, "Crystal Mountain, WA").click()


















