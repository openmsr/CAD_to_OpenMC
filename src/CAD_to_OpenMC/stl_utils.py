"""
a set of utility functions that deal with stl-files
"""
import struct

def merge_stl(dest,src, of='binary'):
  """
  function that takes either a string or a list of strings as input and outputs a single binary stl that
  contains all the facets from the input files and assembles them in the output file \"dest\"
  """
  with open(dest,"wb") as outfp:
    outputbuffer=[]
    outputbh=b''
    ntot=0
    for s in src:
      infp=open(s,"rb")
      header=infp.read(80)
      infp.seek(0)
      if (header[:5].decode('ascii')=='solid'):
        infp.close()
        infp=open(s,'r')
        buffer,n = read_stl_ascii(infp)
      else:
        if s==src[0]:
          outputbh=header
        buffer,n = read_stl_bin(infp)
      infp.close()
      outputbuffer.extend(buffer)
      ntot+=n
  write_stl(dest, outputbuffer,of=of,header=outputbh)

def read_stl_ascii(infp):
  line=infp.readline()
  facets=[]
  import re
  for line in infp:
    line=line.strip()
    if line.startswith('facet normal'):
      facet=[float(x) for x in re.split('\s+',line)[2:] ]
    elif line.startswith('vertex'):
      facet.extend([float(x) for x in re.split('\s+',line.strip())[1:] ])
    elif line.startswith('endfacet'):
      facets.append(facet)
  return facets,len(facets)

def read_stl_bin(infp):
  infp.seek(0)
  header=infp.read(80)
  nfaces=struct.unpack('I',(infp.read(4)))[0]
  buffer=infp.read()
  faces=[f for f in struct.iter_unpack('12fH',buffer)]
  return faces,nfaces

def write_stl(fname,buffer,of='bin',header=b'', name=''):
  if (of.startswith('bin')):
    write_stl_bin(fname,buffer,header)
  else:
    write_stl_ascii(fname,buffer,objectname=name)

def write_stl_ascii(fname,faces,objectname=''):
  outfp=open(fname,'w')
  outfp.write(f'solid {objectname}\n')
  for f in faces:
    outfp.write(f'facet normal {f[0]} {f[1]} {f[2]}\n')
    outfp.write( '  outer loop\n')
    outfp.write(f'    vertex {f[3]} {f[4]} {f[5]}\n')
    outfp.write(f'    vertex {f[6]} {f[7]} {f[8]}\n')
    outfp.write(f'    vertex {f[9]} {f[10]} {f[11]}\n')
    outfp.write( '  endloop\n')
    outfp.write( 'endfacet\n')
  outfp.write(f'endsolid {objectname}\n')
  outfp.close()

def write_stl_bin(fname, faces, header=b'', color=None):
  with open(fname,'wb') as outfp:
    outfp.write(header)
    outfp.write(b'\x00'*(80-len(header)))
    outfp.write(struct.pack('I',len(faces)))
    for facet in faces:
      print(len(facet),*facet)
      outfp.write( struct.pack( '12f',*facet ) )
      #if necessary pad with the attribute byte count (15-bit color) - for now set to 0
      if (len(facet)==12):
        outfp.write(struct.pack('H',0))

if __name__=='__main__':
    import sys
    if len(sys.argv)<3:
      print("Usage: stl_utils.py merge <output.stl> <input0.stl> [<input1.stl>] ....")
      exit(0)
    merge_stl(sys.argv[1],sys.argv[2:])
