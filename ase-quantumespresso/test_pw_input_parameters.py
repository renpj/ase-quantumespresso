from pw_input_parameters import *
## Consistency check
for i in NAMELIST.values():
    for j in i:
        if j not in float_keys+bool_keys+int_keys+char_keys:
            print j

for i in float_keys+bool_keys+int_keys+char_keys:
    lp = True
    for j in NAMELIST.values():
        if i in j:
            lp = False
    if lp:
        print i  
print "Total input parameters: ",len(float_keys+bool_keys+int_keys+char_keys)
print "Float", ":",len(float_keys)
print "Integer", ":",len(int_keys)
print "Character", ":",len(char_keys)
print "Bool", ":", len(bool_keys)
for k,v in NAMELIST.items():
    print "Section",k.upper(),':',len(v)

