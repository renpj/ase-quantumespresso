from pw import *

a = PWInput.from_file('../test_files/pw_scf.in')
a.write_file('test.in')

#a = PWInput.from_file('../test_files/pw_scf_cell.in')
#a.write_file('test_cell.in')
