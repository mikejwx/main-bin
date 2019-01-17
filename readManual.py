#a function that reads the manual classifications
import csv

def readManual(version = 2):
    "reads the manual classifications for 2012."
    "version can be 1, the original (stricter) classification"
    "or version can be 2, the new (less strict) classification"
    
    if version == 1:
        print 'You have requested the older version, this is ill advised'
        manual_CT = []
        manual_NT = []
        manual_OB = []
        with open('/home/xb899100/main/CloudTrailVis_Testing/Manual/training_CT.csv', 'rb') as csvfile:
            spamreader = csv.reader(csvfile, delimiter = ',')
            for row in spamreader:
                manual_CT.append(row)
        with open('/home/xb899100/main/CloudTrailVis_Testing/Manual/training_NT.csv', 'rb') as csvfile:
            spamreader = csv.reader(csvfile, delimiter = ',')
            for row in spamreader:
                manual_NT.append(row)
        with open('/home/xb899100/main/CloudTrailVis_Testing/Manual/training_OB.csv', 'rb') as csvfile:
            spamreader = csv.reader(csvfile, delimiter = ',')
            for row in spamreader:
                manual_OB.append(row)
    elif version == 2:
        manual_CT = []
        manual_NT = []
        manual_OB = []
        with open('/home/xb899100/main/CloudTrailVis_Testing/Manual/Manual_CT.csv', 'rb') as csvfile:
            spamreader = csv.reader(csvfile, delimiter = ',')
            for row in spamreader:
                manual_CT.append(row)
        with open('/home/xb899100/main/CloudTrailVis_Testing/Manual/Manual_NT.csv', 'rb') as csvfile:
            spamreader = csv.reader(csvfile, delimiter = ',')
            for row in spamreader:
                manual_NT.append(row)
        with open('/home/xb899100/main/CloudTrailVis_Testing/Manual/Manual_OB.csv', 'rb') as csvfile:
            spamreader = csv.reader(csvfile, delimiter = ',')
            for row in spamreader:
                manual_OB.append(row)
    else:
        print 'The version you requested does not exist'
    
    return manual_CT[0], manual_NT[0], manual_OB[0]