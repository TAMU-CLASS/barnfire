#! /usr/bin/env python

'''
Andrew Till
Summer 2014

Gets ENDF vii1 XS from T2's site.

For thermal XS, see http://t2.lanl.gov/nis/data/data/ENDFB-VII-thermal/UinUO2' etc.
'''

class Nuclide:
    def __init__(self, name, Z, A):
        self.name = name
        self.Z = Z
        self.A = A

nuclideList = []
nuclideList.append(Nuclide('Pu',94,241))
nuclideList.append(Nuclide('Pu',94,242))
nuclideList.append(Nuclide('U',92,233))

# Print out a list of work to do
with open('../dat/endf/work_items.txt', 'w') as fid:
  for nuclide in nuclideList:
      outname = 'endf_{0}{1}_vii1'.format(nuclide.Z, nuclide.A)
      fid.write('wget http://t2.lanl.gov/nis/data/data/ENDFB-VII.1-neutron/{0}/{1} -O {2}\n'.format(nuclide.name, nuclide.A, outname))

# Print out a script to do the work in parallel
with open('../dat/endf/get_xs.sh', 'w') as fid:
  fid.write('#! /usr/bin/env bash\n')
  fid.write('#wget endf files in parallel!\n')
  fid.write('\n')
  fid.write('cat work_items.txt | parallel\n')
