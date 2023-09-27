#object names from source code
import csv
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By

f = open('AllAntibioticInfo.rtf')
temp = f.readlines()
temp = temp[6::]

antibioticCats = []
antibioticNames = []

for i in temp:
    a = i.find('="')
    b = i.find('">')
    c = i.find('</')
    antibioticCats.append(i[a+2:b])
    antibioticNames.append(i[b+4:c])

f.close()

url = 'http://ecdc.europa.eu/en/healthtopics/antimicrobial_resistance/esac-net-database/Pages/Antimicrobial-consumption-rates-by-country.aspx'

resultid = 'VisibleReportContentctl00_m_g_08b61bd1_a361_4fac_9c75_7c3325705205_ctl00_rv_ctl10'
button = 'ctl00$m$g_08b61bd1_a361_4fac_9c75_7c3325705205$ctl00$btnFilter'

year = 'ctl00$m$g_08b61bd1_a361_4fac_9c75_7c3325705205$ctl00$ddl_Year'
drug = 'ctl00$m$g_08b61bd1_a361_4fac_9c75_7c3325705205$ctl00$ddl_Antimicrobial'

# re-open browser at every loop, in case this helps with connection problems.

filename = 'AllConsumptionData2.csv'

for j in range(70,len(antibioticCats)):
    driver = webdriver.Firefox()
    driver.get(url)
    ofinterest = antibioticCats[j]
    driver.find_element(drug,ofinterest)
    all_data = []
    for i in range(1998,2015):
        driver.find_element(year,str(i))
        # sometimes this fails (not sure why) - handle using exception
        while True:
            try:
                driver.find_element(By.NAME,button).first.click()
                break
            except:
                print ('Problem with clicking, trying again')
        data = [0]
        # for some reason, click sometimes fails - tyr click until manage to get data
        while len(data) < 20:
            data = driver.find_element(By.ID,resultid).text.split('\n')
            print (len(data))
        data = data[4:data.index('United Kingdom')+2]
        chunks=[[antibioticNames[j]]+[str(i)]+data[x:x+2] for x in xrange(0, len(data), 2)]
        all_data.extend(chunks)
    print (j)
    with open(filename, "a") as f:
        writer = csv.writer(f)
        writer.writerows(all_data)
    f.close()
    driver.quit()
exit()