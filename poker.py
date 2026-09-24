from collections import defaultdict
from email.message import EmailMessage
import smtplib
from datetime import date
import os
from dotenv import load_dotenv

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
    payouts = []
    while len(sorted)>0:

        # call the biggest loser "loser", and the biggest winner "winner"

        loser = sorted[0]
        winner = sorted[-1]
        # each of the amounts is divided by a thousand because we multiplied by a thousand earlier
        # case 1: loser loses more than the winner wins
        if abs(loser.amount) > winner.amount:
            loser.amount = loser.amount + winner.amount
            sorted.remove(winner)
            payouts.append(str(loser.name)+' pays '+ str(winner.name)+ ' $' + str(f"{winner.amount/1000:.2f}")+ "\n")

        # case 2: loser loses less than the winner wins    
        elif abs(loser.amount) < winner.amount:
            winner.amount = winner.amount + loser.amount
            sorted.remove(loser)
            payouts.append(str(loser.name)+' pays '+ str(winner.name)+ ' $'+ str(f"{abs(loser.amount)/1000:.2f}")+"\n")

        # case 3: winner and loser are equal    
        elif abs(loser.amount) == abs(winner.amount):
            sorted.remove(winner)
            sorted.remove(loser)
            payouts.append(str(loser.name)+' pays '+ str(winner.name)+ ' $'+ str(f"{winner.amount/1000:.2f}")+"\n")

        # re-sorting before we calculate the next transaction
        unsorted = sorted
        del sorted
        sorted = sort_low_to_high(unsorted)

    print()

    # sorting the payouts so debtors are linked together
    grouped_transactions = defaultdict(list)

    for transaction in payouts:
        name = transaction.split(" pays ")[0]
        grouped_transactions[name].append(transaction)

    debtors = grouped_transactions.keys()

    sorted_payouts = []

    for debtor in debtors:
        sorted_payouts.extend(grouped_transactions[debtor])

    sorted_payouts_onestring = ""
    for items in sorted_payouts:
        sorted_payouts_onestring += items

    # writing outputs to text file
    f1 = open("Ledger.txt", "w")
    f1.writelines(sorted_payouts)

    # setting it up so it sends you an email with the ledger once you run

    load_dotenv()

    SENDER_EMAIL = os.getenv("GMAIL_EMAIL")
    APP_PASSWORD = os.getenv("GMAIL_PW")
    RECIPIENT_EMAIL = input("what email address should I send the ledger to? ")

    msg = EmailMessage()
    msg["Subject"] = "Poker Ledger " +  date.today().isoformat()
    msg["From"] = SENDER_EMAIL
    msg["To"] = RECIPIENT_EMAIL
    msg.set_content("Poker ledger generated on " + date.today().isoformat() + "\n\n" + sorted_payouts_onestring)

    try:
        with smtplib.SMTP("smtp.gmail.com", 587) as smtp:
            smtp.starttls()
            smtp.login(SENDER_EMAIL, APP_PASSWORD)
            smtp.send_message(msg)

            print("email sent successfully")

    except Exception as e:
        print(f"failed to send email. Error: {e}")

    # close the files because it's good practice
    f1.close()
    text.close()

# this else statement runs if the total wasn't equal to 0

else:
    print('inputs do not net to 0. check your numbers and try again')



