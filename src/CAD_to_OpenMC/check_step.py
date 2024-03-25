#a set of function to report on input step-files
import re

def assemble_cmd(line,file):
  line=line.strip()
  while not line.endswith(';'):
    line=line+file.readline().strip()
  return line

def parse_command(command):
  pattern='(#[0-9]+)\s*=\s*(\w+)\s*\((.*)\);'
  multipattern='(#[0-9]+)\s*=\s*\(((\w+)\s*\((.*)\))+\);'
  m=re.match(pattern,command)
  if(not m):
    m=re.match(multipattern,command)
    if (not m):
      return None
  return m.groups()

degen_pattern=r'(#\d+)=DEGENERATE_TOROIDAL_SURFACE\(\'.*\',(#\d+),([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?),([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?),\.([TFtf])\.'

def has_degenerate_toroids(infile='geometry.step', only_known_fail=False):
    dtors=[]
    with open(infile,'r') as f:
        for line in f:
            m=re.search(degen_pattern,line)
            if m:
                dtors.append(m.groups())

    if only_known_fail:
        #remove those toroids which we know work - i.e. those where select outer is true
        for tor in dtors:
            r_maj, r_min = tor[2:4]
            if tor[4] in ['F','f']:
              outer=False
            else:
              outer=True
            if r_maj > 2*r_min or not outer:
                dtors.remove(tor)
    if len(dtors)>0:
        return True, len(dtors)
    else:
        return False, 0
