###Reads the output of CTID_vis.py so that I don't have to run that script
### every time I want to plot new figures.

import csv

def getResults(path = '/home/xb899100/main/CloudTrails/CTID_vis_dates/'):
    testDates = []
    testDates2 = []
    testDates3 = []
    indexes = []
    indexes2 = []
    indexes3 = []
    with open(path + 'new_testDates.csv', 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter = ',')
        for row in spamreader:
            testDates.append(row)
    with open(path + 'new_testDates2.csv', 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter = ',')
        for row in spamreader:
            testDates2.append(row)
    with open(path + 'new_testDates3.csv', 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter = ',')
        for row in spamreader:
            testDates3.append(row)
    with open(path + 'new_indexes.csv', 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter = ',')
        for row in spamreader:
            indexes.append(row)
    with open(path + 'new_indexes2.csv', 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter = ',')
        for row in spamreader:
            indexes2.append(row)
    with open(path + 'new_indexes3.csv', 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter = ',')
        for row in spamreader:
            indexes3.append(row)
    testDates = testDates[0]
    testDates2 = testDates2[0]
    testDates3 = testDates3[0]
    indexes = indexes[0]
    indexes2 = indexes2[0]
    indexes3 = indexes3[0]
    return testDates, testDates2, testDates3, indexes, indexes2, indexes3
