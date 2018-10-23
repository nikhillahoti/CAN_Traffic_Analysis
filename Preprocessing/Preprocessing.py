def preProcessData(filename):

    with open(filename) as f:
        content = f.readlines()

    content = [x.strip() for x in content]

    fIdhData = open("IDH-DATA.txt", "w")
    fIdlData = open("IDL-DATA.txt", "w")
    fIdhIdlData = open("IDH-IDL-DATA.txt", "w")
    fIdhIdl = open("IDH-IDL.txt", "w")
    fData = open("DATA.txt", "w")

    fIdhData.write("IDH,DATA\n")
    fIdlData.write("IDL,DATA\n")
    fIdhIdlData.write("IDH,IDL,DATA\n")
    fIdhIdl.write("IDH,IDL\n")
    fData.write("DATA\n")

    for i in range(len(content)):
        parts = content[i].split(",")
        IDH = parts[0].split(":")[1].strip()
        IDL = parts[1].split(":")[1].strip()
        LEN = parts[2].split(":")[1].strip()
        DATA = parts[3].split(":")[1].replace(" ", "")

        fIdhData.write(IDH + "," + DATA + "\n")
        fIdlData.write(IDL + "," + DATA + "\n")
        fIdhIdlData.write(IDH + "," + IDL + "," + DATA + "\n")
        fIdhIdl.write(IDH + "," + IDL + "\n")
        fData.write(DATA + "" + "\n")

fileName = "drive.dat"
preProcessData(fileName)

