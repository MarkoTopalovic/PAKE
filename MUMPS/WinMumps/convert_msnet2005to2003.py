import os, shutil, sys

# requirements
# ____________

if os.name != 'nt' : sys.exit( "This script is only used on Windows platforms to convert internal dependency project from VisualStudio2005 to VisualStudio2003")

cvstop = os.getenv( 'CVSTOP' )
if not cvstop : sys.exit( "The environment variable CVSTOP has to be set")
if not cvstop in sys.path : sys.path.append( cvstop )
from buildv2.platform import *

for compiler, language in zip( [ F_COMPILER, C_COMPILER, CXX_COMPILER ], [ 'Fortran', 'C', 'C++' ] ) :
  if not compiler.has_key( 'version' ) or not compiler['version']: 
    print dir( compiler )
    sys.exit( "The version of the " + language + "-compiler has to be specified in the platform configuration script" )
if C_COMPILER['version'] != CXX_COMPILER['version'] : sys.exit( "The C and C++ compilers have to be have the same version : " + C_COMPILER['version'] + '|' + CXX_COMPILER['version'] )

# main function
# _____________

def convert() :
  visualstudio = os.path.join( '.', 'VisualStudio' )

  for item in os.listdir( visualstudio + '2005' ) :
    item = os.path.join( visualstudio + '2005', item )
    if os.path.isfile( item ) :
      ext = os.path.splitext( item )[-1]
      if ext == '.vcproj' : vcproj_from_visual_studio_2005_to_2003( item, visualstudio + '2003', float( CXX_COMPILER['version'] ) ) 
      elif ext == '.vfproj' : vfproj_from_visual_studio_2005_to_2003( item, visualstudio + '2003', float( F_COMPILER['version'] ) ) 
      elif ext == '.sln' : solution_from_visual_studio_2005_to_2003( item, visualstudio + '2003' )
 
# visual studio projects
# __________________________

def from_visual_studio_2005_to_2003( filename, output, version ) :
  r = open( filename, 'r' )
  w = open( os.path.join( output, os.path.basename( filename ) ), 'w' )
  stream = ' '.join( r.readlines() )
  
  key = 'Version="'
  if stream.find(key) >= 0 :
    b = stream.find(key)+len(key)
    e = b + stream[b:].find( '"' )
    v = float( stream[b:e].replace(',','.') )
    if v < 8.0 : sys.exit( "The visual studio project " + filename + ' is from VisualStudio2003 (' + v + ')' )
    stream = stream[:b] + str( version ) + stream[e:]
  else : sys.exit( "No version specified in " + filename )
  
  ( key, platform ) = ( '<Platforms>', '' )
  if stream.find(key) >= 0 :
    b = stream.find(key)+len(key)
    e = b + stream[b:].find( '</Platforms>' )
    buffer = repr( stream[b:e] ).replace( '/>', '' ).replace( r'\t', '' ).replace( r'\n', '' ).replace( "'", '' ).replace( ' ', '' ).split( '<Platform' )
    for item in buffer :
      if item.find( 'Win32' ) >= 0 :
	if platform : sys.exit( "Several Win32 platforms are set in " + filename )
        platform = '\t<Platforms>\n\t\t<Platform Name="Win32"/>\n\t</Platforms>'
    if not platform : sys.exit( "No Win32 platform specified in " + filename )
    stream = stream[:b-len(key)-1] + platform + stream[e+len('</Platforms>') :]
  else : sys.exit( "No platform specified in " + filename )


  key = '<Configurations>'
  if stream.find(key) >= 0 :
    b = stream.find(key)+len(key)
    e = b + stream[b:].find( '</Configurations>' )
    buffer = stream[b:e]

    configurations = [ '' ]
    while buffer.find( '<Configuration' ) >= 0 :
      i = buffer.find( '<Configuration' )
      j = buffer.find( '</Configuration>' ) + len( '</Configuration>' )
      config = buffer[i:j]
      buffer = buffer[:i] + buffer[j:]
     
      ( k, l ) = ( config.find( 'Name="' ), len( 'Name="' ) ) 
      if k >= 0 and config[k+l:k+l+config[k+l:].find('"')].find( 'Win32' ) >= 0 : 
        configurations.append( '\t' + config )
    configurations.append( '' )
    stream = stream[:b] + '\n\t'.join( configurations ) + stream[e:]


  key = '<Files>'
  if stream.find(key) >= 0 :
    b = stream.find(key)+len(key)
    e = b + stream[b:].find( '</Files>' )
    buffer = stream[b:e]
    
    k = 0
    while buffer.find( '<FileConfiguration' ) >= 0 :
      i = buffer.find( '<FileConfiguration' )
      j = buffer.find( '</FileConfiguration>' ) + len( '</FileConfiguration>' )
      config = buffer[i:j]
      buffer = buffer[j:]
      u = stream.find( config )
      v = u + len( config )
      if config.find( "Win32" ) < 0 : stream = stream[:u] + stream[v:]
   

    #for k in range( len( buflist ) ) : buflist[k] = '<Filter ' + buflist[k]
    #if r.name.find( 'blacs_mpi_CallFromC' ) >= 0 : print buflist ; sys.exit()


  w.write( stream )
  r.close()
  w.close()
  print '      from ' + r.name + ' to ' + w.name

# C++ visual studio projects
# __________________________

def vcproj_from_visual_studio_2005_to_2003( filename, output, version ) :
  if not os.path.isfile( filename ) and not os.path.splitext( filename )[-1] == '.vcproj' : sys.exit( filename + " is not a C++ visual studio project" )
  else : print '...' + os.path.basename( filename ) + '...'
  from_visual_studio_2005_to_2003( filename, output, version )

# Fortran visual studio projects
# ______________________________

def vfproj_from_visual_studio_2005_to_2003( filename, output, version ) :
  if not os.path.isfile( filename ) and not os.path.splitext( filename )[-1] == '.vfproj' : sys.exit( filename + " is not a Fortran visual studio project" )
  else : print '...' + os.path.basename( filename ) + '...'
  from_visual_studio_2005_to_2003( filename, output, version )

# visual studio solution
# ______________________

def solution_from_visual_studio_2005_to_2003( filename, output ) :
  if not os.path.isfile( filename ) and not os.path.splitext( filename )[-1] == '.sln' : sys.exit( filename + " is not a visual studio solution" )
  else : print '...' + os.path.basename( filename ) + '...'

  r = open( filename, 'r' )
  w = open( os.path.join( output, os.path.basename( filename ) ), 'w' )
  
  for line in r.readlines() :
    if ( line.find( '#' ) == 0 and line.find( '2005' ) > 0 ) or line.find( 'x64' ) >= 0 :
      pass
    else :
      key = 'Format Version 9.00'
      k = line.find( key )
      if k >= 0 : line = line[:k] + 'Format Version 8.00' + line[k+len(key):]
      w.write( line )

  r.close()
  w.close()

# run as a main script
# ____________________

if __name__ == '__main__' :
  convert()

