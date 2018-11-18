
dictOfpackets = {}
globalcounter = 0
def preProcessData(filename, step ,iter):

    global globalcounter
    with open(filename) as f:
        content = f.readlines()

    content = [x.strip() for x in content]
    counter = 0
    i =  200 * iter
    for i in range(i, len(content)):
        parts = content[i].split(",")
        DATA = parts[3].split(":")[1].replace(" ", "")
        for j in range(0, len(DATA), step):
            current = DATA[j: j + step]
            if not current in dictOfpackets:
                dictOfpackets[current] = len(dictOfpackets) + 1
            f2CharacterPacket.write(str(dictOfpackets[current]) + "\n")
            counter += 1
            globalcounter += 1
            if counter > 200: return
    return


step = 2
f2CharacterPacket = open("DATA2", "w")
fileNames = ['autopark.dat', 'drive.dat', 'idle.dat']
iter = 0
while iter < 26:
    for i in range(len(fileNames)):
        preProcessData(fileNames[i], step, iter)
        print("Number of packets written" + str(globalcounter))
        print("Length of Dictionary" + str(len(dictOfpackets)))
    iter += 1
f2CharacterPacket.close()
print("\n Total characters read -> " + str(globalcounter))
print("\n\n\n Awesome !!! File processing ended !!!")
