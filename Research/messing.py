cum_percent_results = []
test_results = [True, False, True, True, True]

for i in range(1):
    print(i)

for i in range(1,len(test_results)+1):
    cum_victory = 0
    for j in range(i):
        if test_results[j] == True:
            cum_victory +=1
    cum_percent_results.append(cum_victory/i*100)

print(cum_percent_results)