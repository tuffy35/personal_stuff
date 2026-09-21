# PEOPLE KEEP GETTING THEMSELVES FOR SOME REASON

import random

# read the names of my CMs and put them into a list

text = open("TeletuffiesNames.txt","r")
NameData = text.readlines()

# stripping of new line characters

for k in range(len(NameData)):
    NameData[k] = NameData[k].strip("\n")

n = len(NameData)
print(n)
buddylist = NameData.copy()
wasntbuddy = 0

# creating a for loop to iterate through each person in the committee at a time (we leave 3 left over because of odd numbers)


for k in range(len(NameData)-3):
    # the person who we are assigning a buddy will be called person. we will just go in order of the list
    person = NameData[k]
    # from this copy, the person who we are assigning a buddy is removed. this way, they don't get assigned to themself
    # note: we have to make sure they haven't already been removed from buddylist
    if person in buddylist == True:
        wasntbuddy = True
        buddylist.remove(person)
    # a random person is selected from this list to be their buddy
    buddyindex = random.randint(0,n-2)
    buddy = buddylist[buddyindex]
    buddylist.remove(buddy)
    # reduce the counter because someone's been removed
    n -= 1
    # add the person who's getting assigned a buddy back to the list so that they can be someone else's buddy
    # note: only add them back if they were there to begin with
    if wasntbuddy == True:
        buddylist.append(person)
    print(person,'\'s secret buddy is', buddy)
needsbuddy = NameData[26:29]
buddylist.reverse()
for j in range(len(needsbuddy)):
    person = needsbuddy[j]
    buddy = buddylist[j]
    print(person,'\'s secret buddy is', buddy)
