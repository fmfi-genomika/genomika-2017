import sys
import re

filename = sys.argv[1]
min_score = int(sys.argv[2])

f = open(filename, "r")

def get_score(s):
    s2 = s[8:]
    s2 = s2[:s2.find(" ")]
    return int(s2)

while True:
    line = f.readline()
    if line == "":
        break
    
    if line[0] == "a":
        line2 = f.readline()
        line3 = f.readline()
        if get_score(line) > min_score:
            print(line.strip())
            print(line2.strip())
            print(line3.strip())
            print("")