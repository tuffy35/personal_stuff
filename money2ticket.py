# This code takes your fundraising total and computes the number of tickets you have in the fundraising lottery.
import math

total = float(input("what is your fundraising total? "))
tickets = 0

if(total>0):
    tickets += math.floor(total/3500)
    tickets += math.floor((total-3500)/500)

    print("you have earned",tickets,"tickets into the IDC lottery")
else:
    print("wouldn't a negative total mean you stole money from THON?")