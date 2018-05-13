# Sample program to test pandas.DataFrame
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

pops = ["Cat","Tiger"]

isin = df['Species'] != 'Cat'
print (isin)
exit()

df1 = df[df['Species'] != 'Cat']
print(df1)

df2 = df[df['Species'].isin(pops)]
print (df2)

# df3 = df
# df3 = df[true/false dataframe of the same length as df]
