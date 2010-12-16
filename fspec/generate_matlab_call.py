#!/usr/bin/env python

import sys
fspec = sys.argv[1]
f = open(fspec+".fspec","r")

data = f.readlines()
data.pop(0)
data.pop(0)

print "lusol_mex('{0}', ...".format(fspec)

for var in data:
    var_name = var.split()[0]
    print "obj." + var_name + ", ..."

## out_str = "lusol_mex('" + fspec + "',"
## for var in data:
##     var_name = var.split()[0]
##     out_str = out_str + "obj." + var_name + ","

## out_str = out_str[:-1]
## out_str = out_str + ");"
## print out_str
