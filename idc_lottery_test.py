import math
import random
import csv
import matplotlib.pyplot as plt
import numpy as np



num_idcs_chosen = 70
num_lotteries = 1000


test_results = []
# each idc is given a number. i populate idc_tickets with that IDC's number as many times as tickets they have. So if IDC 4 has 3 tickets, there are 3 4s in the array.
for i in range(num_lotteries):
    idc_tickets = [] # this array 
    num_tix = 0
    with open("idc_totals_2025.csv",newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for rows in reader:
            idc_number = int(rows["IDC"])
            num_tix = int(rows["Tickets"])
            for i in range(num_tix):
                idc_tickets.append(idc_number)
            
    chosen_idcs = []
    random_index = 0

    while len(chosen_idcs) < num_idcs_chosen + 1:
        random_index = random.randrange(len(idc_tickets))
        idc_it = idc_tickets[random_index]
        chosen_idcs.append(idc_it) # add this idc to the chosen list
        i = 0
        while i < len(idc_tickets):
            if idc_tickets[i] == idc_it:
                idc_tickets.pop(i)   # remove element at index i
                # do NOT increment i here
            else:
                i += 1

    test_flag = False
    for j in chosen_idcs:
        if j == 143:
            test_flag = True
    test_results.append(test_flag)

num_victory = 0 
for k in test_results:
    if k == True:
        num_victory += 1

percent_victory = num_victory/len(test_results)*100

print("You got chosen in",percent_victory,"% of",num_lotteries, "lottery simulations when",num_idcs_chosen, "IDCs were chosen")

cum_percent_results =[]

for i in range(1,len(test_results)+1):
    cum_victory = 0
    for j in range(i):
        if test_results[j] == True:
            cum_victory +=1
    cum_percent_results.append(cum_victory/i*100)

cumresults_plot = np.array(cum_percent_results)
simulations_plot = np.array(range(len(cum_percent_results)))

plt.plot(simulations_plot,cumresults_plot)
plt.xlabel('Lottery Simulation Number')
plt.ylabel('Cumulative Percentage of Winning Lotteries')
plt.title('IDC Lottery Simulation when '+ str(num_idcs_chosen)+ ' IDCs Are Chosen')
plt.show()