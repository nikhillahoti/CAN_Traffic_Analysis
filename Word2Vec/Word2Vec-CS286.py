
from gensim.test.utils import common_texts, get_tmpfile
from gensim.models import Word2Vec

import time
import os

def createDataFiles():
    import time
    start = time.time()

    fileNamesActual = ['autopark.dat', 'drive.dat', 'idle.dat']
    fileOutput = open("Word2VecData/Word2VecDataActual.txt", "w")

    for i in range(len(fileNamesActual)):
        with open(fileNamesActual[i]) as f:
            content = f.readlines()

        # Remove the extra spaces from the sentences
        content = [x.strip() for x in content]
        value = int(0.8 * len(content))

        value = len(content)
        print(value)

        cntCntr = 0
        while cntCntr < value:
            iteration = 0
            currMessagePackets = ""

            while iteration < 50 and cntCntr < value:
                parts = content[cntCntr].split(",")

                # this is the whole message
                DATA = parts[3].split(":")[1].replace(" ", "")

                currMessagePackets += str(DATA) + " "

                cntCntr += 1
                iteration += 1

            if cntCntr < value:
                cntCntr -= 40
            else:
                break
            fileOutput.write(currMessagePackets + "\n")
    fileOutput.close()

    fileOutput = open("Word2VecData/Word2VecDataSimulated.txt", "w")
    fileNamesSimulated = ['New_DRIVE_Data.txt', 'New_IDLE_Data.txt', 'DRIVE.rtf', 'IDLE.rtf']

    for i in range(len(fileNamesSimulated)):
        with open(fileNamesSimulated[i]) as f:
            content = f.readlines()

        # Remove the extra spaces from the sentences
        content = [x.strip() for x in content]
        print(len(content))

        cntCntr = 0
        while cntCntr < len(content):
            iteration = 0
            currMessagePackets = ""

            while iteration < 50 and cntCntr < len(content):
                DATA = content[cntCntr][15:38]

                # this is the whole message
                parts = DATA.split()

                tempo = ""
                for j in range(8):
                    if j >= len(parts) or parts[j] == "  ":
                        tempo += "00"
                    else:
                        tempo += parts[j]

                currMessagePackets += str(tempo) + " "
                cntCntr += 1
                iteration += 1

            if cntCntr < len(content):
                cntCntr -= 40
            else:
                break
            fileOutput.write(currMessagePackets + "\n")

    fileOutput.close()
    end = time.time()

    print("Awesome !!! File processing done !!!")
    print("Total Time for file processing ---> ", end - start)

class IteratingClass:
    def __init__(self, dirName):
        self.dirName = dirName

    def __iter__(self):
        for fName in os.listdir(self.dirName):
            for line in open(os.path.join(self.dirName, fName)):
                yield line.split()

def createWord2VecModels():
    start = time.time()
    dataDirec = IteratingClass('/home/nikhil/PycharmProjects/CodingPrac/Word2VecData')
    model = Word2Vec(dataDirec, size=200, window=5, min_count=1, workers=8)
    model.save(fileName)
    end = time.time()
    print("\n\nTraining Successful for Word2Vec Model!!!")
    print("Total Time for Word2Vec model -> ", (end - start))

def createDNNModel():
    pass


def train():
    createDataFiles()
    createWord2VecModels()

    createDNNModel()

def evaluate():
    model = Word2Vec.load(fileName)
    print("Loaded successfully")
    #print(model.wv['FEB0FF999E500000'])
    print(len(model.wv['FEB0FF999E500000']))

    #print(model.wv['000200000000002A'])
    print(len(model.wv['000200000000002A']))


fileName = "CS286-Word2Vec.model"

train()
#evaluate()

