import random

# read the names of my CMs and put them into a list

text = open("TeletuffiesNames.txt","r")
NameData = text.readlines()

# stripping of new line characters

for k in range(len(NameData)):
    NameData[k] = NameData[k].strip("\n")

# n is the number of people still in the list that haven't been assigned a buddy yet
    
n = len(NameData)

# initializing a groupedlist that will contain sublists of all the pairings

GroupedList = []

# 2 random indices will be selected to be paired together, and then those names will be removed from the list. Then, the range of
# random indices will shorten to the length of the new list. 

while n > 1:
    index1 = random.randint(0,n-1)
    person1 = NameData.pop(index1)
    index2 = random.randint(0,n-2)
    person2 = NameData.pop(index2)
    NewGroup = [person1, person2]
    GroupedList.append(NewGroup)
    n = n-2

# because there is an odd number, we will throw the last remaining person in with the last group in the list
    
lastperson = NameData.pop(0)
GroupedList[-1].append(lastperson)
print(GroupedList)

text.close()