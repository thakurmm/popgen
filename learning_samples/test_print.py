# The equivalent of this in python 2.7
# print ("arg %d") % 10
# is as below
anystr = "arg {0}".format(10)
print(anystr)

# int is data type 'd'
int_anystr = "arg {0:d}".format(10)
print (int_anystr)

# float is data type 'f'
float_anystr = "arg {0:f}".format(10)
print (float_anystr)
