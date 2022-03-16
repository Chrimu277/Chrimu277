print("Hello World")
a = 1
b = 7
c = "asd"
e = a+b
print(e)
# this is a comment

import numpy as np


#ARRAYS (lists)
print("\n-----ARRAYS-----")
array1 = [1,2,7,6,-1,13]
friends = ["Kevin","Jim","Karen","Jim"]
print(array1)
print(friends)
print(friends[0])
friends.extend(["JFK","OBAMA","McCain"])
friends.insert(1,"Chris")
friends.append(array1)
print(friends)
friends.remove("JFK")
friends.pop()
print(friends.index("McCain"))
print(friends.count("Jim"))
friends.sort()
array1.sort()
print(friends)
print(array1)

#TUPLE, (immutable list)
print("\n-----TUPLES-----")
coordinates = (2, 5)
coord2 = 2,5
print(coord2)
print(coordinates)
print(coordinates[0])


#LOOPS
print("\n-----LOOPS-----")

#class


firstname = "John"
lastname = "Farse"
isMale = True
age = 35

fullname = firstname+" "+lastname


#PRINTING
print("\n-----PRINTING-----")
print("Once there was a man named "+fullname)
print("I am "+str(age)+" years old")


#STRINGS
print("\n-----STRINGS-----")
print("\t Farse is "+str(len(lastname))+ " letters long")
print(firstname[0] + "."+lastname)
print(lastname.index("se"))
print(lastname.index("F"))
print(lastname.replace("se","miller"))

#NUMBERS
print("\n-----NUMBERS-----")
print(abs(-37.4))
print(pow(2,8))
print(str(max(4,6))+" and "+str(min(2,3)))
print("round 2.7 to "+ str(round(2.7)))

#MATH IMPORT (math is a module, * means all)
print("\n-----MATH MODULE-----")
from math import *
print("-3.2 floors to "+str(floor(-3.2)))
print(sqrt(78))

#USER INPUT
print("\n-----USER INPUT-----")
#username = input("Please enter your name:\n--> ")
#userage = input("Please enter your age:\n--> ")
#print("Hello, "+username+". You apparently are "+str(userage)+" years old")

#FUNCTIONS
print("\n-----FUNCTIONS-----")
def myfunction():
    print("First function")

def cubeit(numb):
    return numb*numb*numb

print(cubeit(4))
myfunction()

#IF STATEMENT
print("\n-----IF STATEMENT-----")
if isMale:
    print("he is male")
else:
    print("she is female")

isTall = True

if isMale or isTall:
    print("You must be strong arent you")
else:
    print("You are wether tall nor male")

if isMale and isTall:
    print("You are even male and strong !!")

#CASTS
print("\n-----CASTS-----")
print(int(3.2))
print(str(-4.3))
print(float("3,2"))
print(float("3.2"))
print(float("-1"))