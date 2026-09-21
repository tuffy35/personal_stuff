# importing data from text file

data = open("NameList.txt","r")

# separating data line by line

data = data.readlines()

# removes \n, separates each line at the comma, and converts all numbers from strings into floats

def process(data):

    print(data)
    for k in range(len(data)):
        data[k] = data[k].strip()
        data[k] = data[k].split(", ")
        data[k][1] = float(data[k][1])
    print(data)

data2 = ['Nock, 30\n', 'Moke, 35\n', 'Androw, 4']
process(data2)


