import math

tickets = int(input("How many tickets do you want? "))
total = 0

if tickets > 0:
    # one ticket for 3500
    total += 3500
    tickets -=1

    count = 0 # this is how we'll determine when to subtract an extra ticket for each 3500
    while tickets>0:
        total+=500
        tickets-=1
        count +=1
        if count == 7:
            tickets-= 1
            count = 0

    print("you need to raise at least $",total)

else:
    print("now why would you want",tickets,"tickets")