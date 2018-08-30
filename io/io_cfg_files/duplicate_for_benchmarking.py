#!/usr/bin/python

import sys

num_data=2

def process_diagnostics(diag_lines, new_file):
    field_name=diag_lines[0].split("field=\"")[1].split("\"")[0].split("_")[0]
    field_name_sc=field_name.swapcase()
    for i in range(1,num_data):
        for line in diag_lines:
            updated=line.replace(field_name, field_name+"~"+str(i))
            updated=updated.replace(field_name_sc, field_name_sc+"~"+str(i))
            new_file.write(updated)

inside_dd=False
inside_diagnostics=False

diag_list=[]

new_file = open(sys.argv[1]+"_updated","w") 

with open(sys.argv[1]) as fp:  
   line = fp.readline()   
   while line:
       new_file.write(line)
       if "<data-definition" in line:
           inside_dd=True
       if "</data-definition>" in line:
           inside_dd=False
       if "<diagnostic" in line:
           inside_diagnostics=True
       if (inside_diagnostics):
           diag_list.append(line)
       if "</diagnostic>" in line:
           inside_diagnostics=False
           process_diagnostics(diag_list, new_file)
           del diag_list[:]       
       if ("<member" in line):
           field_name=line.split("name=\"")[1].split("\"")[0].split("_")[0]
           for i in range(1,num_data):
               updated=line.replace(field_name, field_name+"~"+str(i))
               new_file.write(updated)
       if (inside_dd and "<field name" in line):
           name=line.split()[1].split("\"")[1]
           first_quote=line.find("\"")+1
	   restofline=line[first_quote+line[first_quote:].find("\"")+1:]
           for i in range(1,num_data):
              new_file.write("<field name=\"" + name+"~"+str(i)+"\""+restofline)
       line = fp.readline()

new_file.close()