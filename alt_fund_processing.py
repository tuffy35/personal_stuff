#######################################################################################
#   This is a short python script that I used to calcualte the grand totals for       #
#   the THON 2024 alternative fundraising data. In order for this to work             #
#   properly, you'll need to have the alt fund data in a .csv called "alt_funds_2024" #
#   as well as a blank .csv called "alt_funds_2024_totals", and they must both be in  #
#   the same directory as this file. Please don't roast me if this is inefficient or  #
#   something, I learned python from a textbook and a few google searches. If you     #
#   want to iterate on it feel free! Just keep a log of the version and the edits     #
#   made.                                                                             #
#######################################################################################
#   Nick Tufillaro, FOC  #
#   7/3/2024             #
#   V1.0                 #
##########################



import csv
from collections import defaultdict

# initializing a default dictionary with keys = org type, values = total amount credited
fundraising_totals = defaultdict(int)

# opening the alt funds datasheet
with open("alt_funds_2024.csv",newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for rows in reader:
        if rows["credited"] == '':  # if no amount was credited, we skip that fundraiser
            continue
        else:
            credited = int(float(rows["credited"])*100)     # converting the credited amount to an integer and multiplying by 100 to avoid float error 
        org_type = rows["org_type"]
        fundraising_totals[org_type] += credited    # adding credited amount to the org type counter

csvfile.close()

# formatting the values such that they are to exactly 2 decimal places
for key, value in fundraising_totals.items():
    fundraising_totals[key] = format(value/100,"0.2f")

# calculating the grand total so we can add that as a row
total = 0
for items in fundraising_totals.items():
    total += float(items[1])
total = format(total,"0.2f")

# writing the values to a new csv file
with open("alt_funds_2024_totals.csv",'w',newline='') as csvfile:
    fieldnames = ["org type","fundraising total"]
    writer = csv.writer(csvfile)
    writer.writerow(fieldnames)     # writing out the headers
    for items in fundraising_totals.items():    # writing out each row from the dictionary values
        writer.writerow([items[0],items[1]])    # .items() method converts each dictionary key-value pair into a tuple, so we can index with a 0 and 1
    writer.writerow(["total",total])  # writing the total

csvfile.close()
