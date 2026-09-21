# NEXT UPDATE: after payments are all done, organize them such that each debtor's
# payments are listed consecutively to improve readability
# creating a person class that contains a person's name and the amount they are up/down

class Person:
    def __init__(self, name, amount):
        self.name = name
        self.amount = amount

# creating function to sort the arrays
def sort_low_to_high(unsorted):
    sorted = []
    while len(unsorted)>0:
        SmallIndex = 0
        for m in range(len(unsorted)):
            if unsorted[m].amount < unsorted[SmallIndex].amount:
                SmallIndex = m
        sorted.append(unsorted.pop(SmallIndex))
    return sorted

# uploading and reading the data from a text file named "NameList"     

text = open("NameList.txt","r")
data = text.readlines()

# reformatting data by stripping of new line characters and spaces, separating names and numbers, and converting amounts to floats

for k in range(len(data)):
    data[k] = data[k].strip()
    data[k] = data[k].split(", ")
    # multiplying the amount by 1000 and then converting to int to prevent float error
    data[k][1] = int(float(data[k][1])*1000)

# check if the numbers add up to 0 before running the entire program

numbers = [i[1] for i in data]
total = 0
for n in numbers:
    total += n
print('check: total =', total/1000.,',')
if total == 0:
    print('so the values are good \n')
    # organizing names and amounts into Person class objects in an unsorted list. Each person is assigned a name and an amount

    unsorted = [*range(len(data))]
    for j in range(len(data)):
        unsorted[j] = Person(data[j][0], data[j][1])

    # if anyone's net winnings is 0, remove them from the list because they don't need to pay anyone or get paid anything

    for p in range(len(unsorted)):
        if unsorted[p].amount == 0:
            del unsorted[p]

    # initial sort
    sorted = sort_low_to_high(unsorted)

    # loop through the sorted list, making the biggest loser pay the biggest winner. continue to work in until no more money is owed

    while len(sorted)>0:

        # call the biggest loser "loser", and the biggest winner "winner"

        loser = sorted[0]
        winner = sorted[-1]
        # each of the amounts is divided by a thousand because we multiplied by a thousand earlier
        # case 1: loser loses more than the winner wins
        if abs(loser.amount) > winner.amount:
            loser.amount = loser.amount + winner.amount
            sorted.remove(winner)
            print(loser.name,'pays', winner.name, '$', f"{winner.amount/1000:.2f}")

        # case 2: loser loses less than the winner wins    
        elif abs(loser.amount) < winner.amount:
            winner.amount = winner.amount + loser.amount
            sorted.remove(loser)
            print(loser.name,'pays', winner.name, '$', f"{abs(loser.amount)/1000:.2f}")

        # case 3: winner and loser are equal    
        elif abs(loser.amount) == abs(winner.amount):
            sorted.remove(winner)
            sorted.remove(loser)
            print(loser.name,'pays', winner.name, '$', f"{winner.amount/1000:.2f}")

        # re-sorting before we calculate the next transaction
        unsorted = sorted
        del sorted
        sorted = sort_low_to_high(unsorted)

    print()

# this else statement runs if the total wasn't equal to 0

else:
    print('your doodoo brain inputted the numbers wrong')

# close the file because it's good practice

text.close()
