""" M3C 2016 Homework 1
Python script which prints name and college id
To run this code, enter "python hw1.py" at the terminal
"""

#1. modify the list, Output, so that it contains your name and college id
Output = ["Evgeniia Gleizer","00948999"]

#2. modify x and y in the print statements below so that your name and college id are output
print "M3C 2015 Homework 1 by", Output[0]
print "CID:", Output[1]


#3. Add python code here which removes any leading zeros from your CID and stores the result in
#   the variable, ID2. For example, "00000001" would become "1"
ID2=Output[1]
while ID2[0]=='0':
    ID2=ID2[1:len(ID2)]
print "ID2:", ID2

#4. Add python code here which removes all zeros from your CID and stores the result in
#   the string, ID3. For example, "10000001" would become "11"
ID3=str()
for x in Output[1]:
    if x>"0":
        ID3+=str(x)
        
print "ID3:", ID3


#Note: Your code for parts 3 and 4 should work for any 8-digit CID with at least one non-zero
#number. 
