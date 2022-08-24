#!/usr/bin/env python
import sys
import re
import os
sys.setrecursionlimit(100000)

argv=sys.argv

scaling_factor=1
if ('--m2in') in argv:
  #data in the step-file is in meters but will be imported into an inches drawing
  scaling_factor=100/2.54
  argv.remove('--m2in')
if('--in2m') in argv:
  scaling_factor=2.54/100
  argv.remove('--in2m')

try:
  stepfile=argv[1]
except:
  print("ERROR: please supply a step file")
  exit(1)

try:
  startlabel=argv[2]
except:
  print('WARNING: no startlabel found - using #1234')
  startlabel="#1234"

def assemble_cmd(line,file):
  line=line.strip()
  while not line.endswith(';'):
    line=line+file.readline().strip()
  return line

def parse_command(command):
  secpattern='([A-Z_\-0-9]+);'
  pattern='(#[0-9]+)\s*=\s*(\w+)\s*\((.*)\);'
  multipattern='(#[0-9]+)\s*=\s*\((:?(\w+)\s*\((.*)\)\s)+\);'
  m=re.match(pattern,command)
  if(not m):
    m=re.match(multipattern,command)
    if m:
      ret=[m.group(1),[],[]]
      m2=re.findall('\w+\(.*\)')
      if m2:
        for l in m2:
          m3=re.match('(\w+)\((.*)\)',l)
          ret[1].append(m3.group(1))
          ret[2].appemd(m3.group(2))
      print(ret)
      return ret
    if (not m):
      m=re.match(secpattern,command)
      if m:
        return [m.groups(1),'','']
      return None
  return m.groups()

def find_dependents_loop(stepfile, label, picked=[]):
  alldeps=[]
  deps=[label]
  done=False
  index=0
  while not done:
    done=True
    index_update=0
    loopdeps=[]
    with open(stepfile) as f:
      for line in f:
        stepcmd=assemble_cmd(line,f)
        splitcmd=parse_command(stepcmd)
        if not splitcmd:
          continue
        if (splitcmd[0] in deps[index:]):
          ls=list(set(re.findall('#\d+',"".join(splitcmd[-1]))))
          index_update=index_update+1
          loopdeps.extend(ls)
          done=False
    #we have now extended the list of labels to be checked for dependencies
    deps.extend(loopdeps)
    index=index+index_update
  return deps

def rfind_dependents_loop(stepfile, label, picked=[]):
  alldeps=[]
  deps=[label]
  done=False
  index=0
  while not done:
    done=True
    loopdeps=[]
    with open(stepfile) as f:
      for line in f:
        stepcmd=assemble_cmd(line,f)
        splitcmd=parse_command(stepcmd)
        if not splitcmd:
          continue
        if (deps[-1] in re.findall('#\d+',"".join(splitcmd[-1]))):
          loopdeps.extend([splitcmd[0]])
          done=False
      #we have now extended the list of labels to be checked for dependencies
    deps.extend(loopdeps)
  return deps

def find_dependents(stepfile,label,picked=[]):
  alldeps=[]
  deps=[]
  with open(stepfile) as f:
    for line in f:
      stepcmd=assemble_cmd(line,f)
      splitcmd=parse_command(stepcmd)
      if not splitcmd:
        continue
      if (label in splitcmd[0:]):
        #split the command argument to print(splitcmd)
        deps.extend(re.findall('#\d+',"".join(splitcmd[-1])))
  #print(deps)
  for d in deps:
    dd=set(d) ^ set(picked)
    alldeps.append(list(dd))
    alldeps.extend(find_dependents(stepfile,list(dd),picked=picked+list(dd)))
    #print(alldeps)
  return alldeps

def rfind_dependents(stepfile,label):
  alldeps=[]
  deps=[]
  with open(stepfile) as f:
    for line in f:
      stepcmd=assemble_cmd(line,f)
      splitcmd=parse_command(stepcmd)
      if not splitcmd:
        continue
      if (label in [l for l in re.split('[(),]',splitcmd[-1]) if re.match('#\d+',l)]):
        #split the command argument to print(splitcmd)
        deps.extend([splitcmd[0]])
  for d in deps:
    alldeps.append(d)
    alldeps.extend(rfind_dependents(stepfile,d))
  return alldeps

def find_line(stepfiles,label=None,cmd='CLOSED_SHELL',arg=None):
  with open(stepfile) as f:
    for line in f:
      foundlabel=foudncmd=foundarg=False
      stepcmd=assemble_cmd(line,f)
      splitcmd=parse_command(stepcmd)
      if label is None or label==splitcmd[0]:
        foundlabel=True
      if label is None or cmd==splitcmd[1]:
        foundcmd=True
      if label is None or arg in splitcmd[2]:
        foundarg=True
      if foundlabel and foundcmd and foundarg:
        return splitcmd

def print_splitcmd(splitcmd=['#1','CARTESIAN_POINT',[0,0,0]]):
  print(f'{splitcmd[0]}={splitcmd[1]}({",".join(splitmcd[2])});')

def ismultiline(cmd):
  m=re.match('#\d+=\(',cmd)
  if m:
    return True
  else:
    return False

def splitmultiline(cmd):
  m=re.match(r'(#\d+=\()((?:\w+\(.+\))+)(\);)',cmd)
  if m:
    start=m.group(1)
    end=m.group(3)
    return  [start]+ re.findall(r'\w+\(.*?\)?\)',m.group(2))+[end]
  else:
    return None

def findprint_lines(stepfiles,labels):
  with open(stepfile) as f:
    for line in f:
      stepcmd=assemble_cmd(line,f)
      m=re.match('(#\d+).*',stepcmd)
      if m is not None:
        s=m.group(1)
        #print(s)
        if s in labels:
          if ismultiline(stepcmd):
            #this is a multiline command should be written as such
            s=splitmultiline(stepcmd)
            #m3=re.match(r'(#\d+=\()((?:\w+\(.+\))+)(\);)',stepcmd)
            #s=m3.group(2)
            #print(m3.group(1))
            print(("\n".join(s)).strip())
            #print(m3.group(3))
          else:
            #single line command
            print(m.group(0))
    return True
  return False

header="""
ISO-10303-21;
HEADER;
/* Generated by software containing ST-Developer
 * from STEP Tools, Inc. (www.steptools.com)
 */
/* OPTION: using custom renumber hook */

FILE_DESCRIPTION(
/* description */ ('STEP AP242',
'CAx-IF Rec.Pracs.---Representation and Presentation of Product Manufa
cturing Information (PMI)---4.0---2014-10-13',
'CAx-IF Rec.Pracs.---3D Tessellated Geometry---0.4---2014-09-14','2;1'),

/* implementation_level */ '2;1');

FILE_NAME(
/* name */ '62f6081ce0a8eb784fad3736',
/* time_stamp */ '2022-08-12T07:58:40+00:00',
/* author */ (''),
/* organization */ (''),
/* preprocessor_version */ 'ST-DEVELOPER v18.102',
/* originating_system */ '  ',
/* authorisation */ '  ');

FILE_SCHEMA (('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));
ENDSEC;

DATA;
"""
footer="""
#395963=(
LENGTH_UNIT()
NAMED_UNIT(*)
SI_UNIT($,.METRE.)
);
#395964=PRODUCT_DEFINITION_SHAPE('','',#395965);
#395965=PRODUCT_DEFINITION('','',#395967,#395966);
#395966=PRODUCT_DEFINITION_CONTEXT('',#395973,'design');
#395967=PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE('','',#395969,
 .NOT_KNOWN.);
#395968=PRODUCT_RELATED_PRODUCT_CATEGORY('','',(#395969));
#395969=PRODUCT('graphite_allstringersin1','graphite_allstringersin1',
'graphite_allstringersin1',(#395971));
#395970=PRODUCT_CATEGORY('','');
#395971=PRODUCT_CONTEXT('',#395973,'mechanical');
#395972=APPLICATION_PROTOCOL_DEFINITION('international standard',
'ap242_managed_model_based_3d_engineering',2011,#395973);
#395973=APPLICATION_CONTEXT('managed model based 3d engineering');
ENDSEC;
END-ISO-10303-21;
"""

def labelkey(label):
  m=re.match(r'#(\d+)',label)
  if m:
    return int(m.group(1))
  else:
    return 10000000


if __name__ == '__main__':
  #initial search for up and down depedencies from the starting label
  up=rfind_dependents_loop(stepfile,startlabel)
  up=list(set(up))
  up.sort(key=labelkey)
  print(up)
  down=[]
  for label in up:
    down1=find_dependents_loop(stepfile,label)
    print(label,"      ",down1)
    down.extend(down1)
  down=list(set(down))
  full=(up+down)
  full=list(set(full))
  full=sorted(full,key=labelkey)
  #print(len(full))
  #we should now have a complete coherent set we can now assembler
  print(header)
  findprint_lines(stepfile,full)
  print(footer)
