# https://www.python-course.eu/python3_lambda.php

# lambda is a python way of defining a function on a single line.
# In the example below, we are defining a function "sum" that takes two arguments and performs the addition operation on them
sum = lambda x,y: x + y
# We can then use this function as below
print (sum(2,3))

# The above operation could also have been done as
# First define a function
def f_sum(x,y):
    return (x+y)

print (f_sum(2,3))

# ------------
# This is used in pandas dataframe in the DataFrame.apply(), by passing it as a function.
import pandas as pd
pet_names = [
    ["Jack", "Cat"],
    ["Jill", "Dog"],
    ["Tom", "Cat"],
    ["Harry", "Dog"],
    ["Hannah", "Dog"],
    ["John", "Tiger"],
    ["Dan", "Tiger"]
]
df = pd.DataFrame(pet_names, columns=["Name", "Species"])

new_df = df.apply(lambda str:str + "_")
print (new_df)
