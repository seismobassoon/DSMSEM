#!/usr/bin/env python

# ./mt_chunk --lr 50 --ur 100 --dr 10 --llat -45 --ulat 45 --dlat 10 --llon -45 --ulon 45 --dlon 10

import sys
import numpy
import math
import h5py

# Usage message
if len(sys.argv) == 1:
    print "usage: ./mt_chunk --lr lr --ur ur --dr dr --llat llat --ulat ulat --dlat dlat --llon llon --ulon ulon --dlon dlon"
    print "mt_chunk is a utility designed to create a shell (HEXA8)"
    print "mandatory arguments are described here after:"
    print "--lr: lower radius"
    print "--ur: upper radius"
    print "--dr: radius step"
    print "--llat: lower latitude in degrees"
    print "--ulat: upper latitude in degrees"
    print "--dlat: distance step along latitude"
    print "--llon: lower longitude in degrees"
    print "--ulon: upper longitude in degrees"
    print "--dlon: distance step along longitude"
    print "optional arguments are described here after:"
    print "--h27: generate HEXA27"
    sys.exit(0)

# Check command line
lr = "none"
ur = "none"
dr = "none"
llat = "none"
ulat = "none"
dlat = "none"
llon = "none"
ulon = "none"
dlon = "none"
h8 = True
for a in range(1, len(sys.argv)):
    arg = sys.argv[a]
    if arg == "--lr":
        lr = float(sys.argv[a+1]) if a+1 < len(sys.argv) else "none"
    if arg == "--ur":
        ur = float(sys.argv[a+1]) if a+1 < len(sys.argv) else "none"
    if arg == "--dr":
        dr = float(sys.argv[a+1]) if a+1 < len(sys.argv) else "none"
    if arg == "--llat":
        llat = float(sys.argv[a+1]) if a+1 < len(sys.argv) else "none"
    if arg == "--ulat":
        ulat = float(sys.argv[a+1]) if a+1 < len(sys.argv) else "none"
    if arg == "--dlat":
        dlat = float(sys.argv[a+1]) if a+1 < len(sys.argv) else "none"
    if arg == "--llon":
        llon = float(sys.argv[a+1]) if a+1 < len(sys.argv) else "none"
    if arg == "--ulon":
        ulon = float(sys.argv[a+1]) if a+1 < len(sys.argv) else "none"
    if arg == "--dlon":
        dlon = float(sys.argv[a+1]) if a+1 < len(sys.argv) else "none"
    if arg == "--h27":
        h8 = False
if lr == "none":
    print "error : missing lr option at command line"
    sys.exit(1)
if ur == "none":
    print "error : missing ur option at command line"
    sys.exit(1)
if lr >= ur:
    print "error : lr >= ur"
    sys.exit(1)
if dr == "none":
    print "error : missing dr option at command line"
    sys.exit(1)
if llat == "none":
    print "error : missing llat option at command line"
    sys.exit(1)
if ulat == "none":
    print "error : missing ulat option at command line"
    sys.exit(1)
if llat >= ulat:
    print "error : llat >= ulat"
    sys.exit(1)
if dlat == "none":
    print "error : missing dlat option at command line"
    sys.exit(1)
if llon == "none":
    print "error : missing llon option at command line"
    sys.exit(1)
if ulon == "none":
    print "error : missing ulon option at command line"
    sys.exit(1)
if llon >= ulon:
    print "error : llon >= ulon"
    sys.exit(1)
if dlon == "none":
    print "error : missing dlon option at command line"
    sys.exit(1)

# Create nodes
xAxis = numpy.arange(llat, ulat+dlat, dlat) if h8 else numpy.arange(llat, ulat+dlat/2., dlat/2.)
yAxis = numpy.arange(llon, ulon+dlon, dlon) if h8 else numpy.arange(llon, ulon+dlon/2., dlon/2.)
rAxis = numpy.arange(lr,   ur  +dr,   dr  ) if h8 else numpy.arange(lr,   ur  +  dr/2.,   dr/2.)
nID = 0
nodes = []
for ix, x in enumerate(xAxis):
    nodes.append([])
    for iy, y in enumerate(yAxis):
        nodes[ix].append([])
        for ir, r in enumerate(rAxis):
            xx = r*math.sin((90-x)*math.pi/180)*math.cos(y*math.pi/180)
            yy = r*math.sin((90-x)*math.pi/180)*math.sin(y*math.pi/180)
            rr = r*math.cos((90-x)*math.pi/180)
            nodes[ix][iy].append({"x" : xx, "y" : yy, "r" : rr, "nID" : nID})
            nID = nID+1

# Create elements
step = 1 if h8 else 2
elems = []
for ix in range(0, len(xAxis)-1, step):
    for iy in range(0, len(yAxis)-1, step):
        for ir in range(0, len(rAxis)-1, step):
            n00 = nodes[ix     ][iy     ][ir     ]["nID"]
            n01 = nodes[ix+step][iy     ][ir     ]["nID"]
            n02 = nodes[ix+step][iy+step][ir     ]["nID"]
            n03 = nodes[ix     ][iy+step][ir     ]["nID"]
            n04 = nodes[ix     ][iy     ][ir+step]["nID"]
            n05 = nodes[ix+step][iy     ][ir+step]["nID"]
            n06 = nodes[ix+step][iy+step][ir+step]["nID"]
            n07 = nodes[ix     ][iy+step][ir+step]["nID"]
            if h8:
                elems.append({"n00" : n00, "n01" : n01, "n02" : n02, "n03" : n03, "n04" : n04, "n05" : n05, "n06" : n06, "n07" : n07})
            else:
                n08 = nodes[ix+1][iy  ][ir  ]["nID"]
                n09 = nodes[ix+2][iy+1][ir  ]["nID"]
                n10 = nodes[ix+1][iy+2][ir  ]["nID"]
                n11 = nodes[ix  ][iy+1][ir  ]["nID"]
                n12 = nodes[ix+1][iy  ][ir+2]["nID"]
                n13 = nodes[ix+2][iy+1][ir+2]["nID"]
                n14 = nodes[ix+1][iy+2][ir+2]["nID"]
                n15 = nodes[ix  ][iy+1][ir+2]["nID"]
                n16 = nodes[ix  ][iy  ][ir+1]["nID"]
                n17 = nodes[ix+2][iy  ][ir+1]["nID"]
                n18 = nodes[ix+2][iy+2][ir+1]["nID"]
                n19 = nodes[ix  ][iy+2][ir+1]["nID"]
                n20 = nodes[ix  ][iy+1][ir+1]["nID"]
                n21 = nodes[ix+2][iy+1][ir+1]["nID"]
                n22 = nodes[ix+1][iy  ][ir+1]["nID"]
                n23 = nodes[ix+1][iy+2][ir+1]["nID"]
                n24 = nodes[ix+1][iy+1][ir  ]["nID"]
                n25 = nodes[ix+1][iy+1][ir+2]["nID"]
                n26 = nodes[ix+1][iy+1][ir+1]["nID"]
                elems.append({"n00" : n00, "n01" : n01, "n02" : n02, "n03" : n03, "n04" : n04, "n05" : n05, "n06" : n06, "n07" : n07,
                              "n08" : n08, "n09" : n09, "n10" : n10, "n11" : n11, "n12" : n12, "n13" : n13, "n14" : n14, "n15" : n15,
                              "n16" : n16, "n17" : n17, "n18" : n18, "n19" : n19, "n20" : n20, "n21" : n21, "n22" : n22, "n23" : n23,
                              "n24" : n24, "n25" : n25, "n26" : n26})

# Write HDF5 file
f = h5py.File("mt_chunk.h5", "w")
d = [(nodes[x][y][r]["x"], nodes[x][y][r]["y"], nodes[x][y][r]["r"]) for x in range(len(xAxis)) for y in range(len(yAxis)) for r in range(len(rAxis))]
dn = f.create_dataset("/Nodes", data=d, dtype=float)
dn.attrs["nn"] = nID
eLabel = "Hexa8" if h8 else "Hexa27"
nPerElem = 8 if h8 else 27
d = [[elems[e]["n%02d"%n] for n in range(nPerElem)] for e in range(len(elems))]
de = f.create_dataset("/Sem3D/"+eLabel, data=d, dtype=int)
de.attrs["ne"] = len(elems)
d = [0 for e in range(len(elems))]
f.create_dataset("/Sem3D/Mat", data=d, dtype=int)

# Write XDMF file
xmf = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">
  <Domain>
    <Grid Name="%(name)s">
      <Geometry GeometryType="XYZ">
        <DataItem Format="HDF" DataType="Float" Rank="2" Dimensions="%(nn)d 3">
          %(name)s.h5:/Nodes
        </DataItem>
      </Geometry>
      <Topology TopologyType="%(eType)s" NumberOfElements="%(ne)d">
        <DataItem Format="HDF" DataType="Int" Rank="2" Dimensions="%(ne)d %(nPerElem)d">
          %(name)s.h5:/Sem3D/%(eLabel)s
        </DataItem>
      </Topology>
    </Grid>
  </Domain>
</Xdmf>
"""
eType = "Hexahedron" if h8 else "Hexahedron_27"
open("mt_chunk.xmf", "w").write(xmf%{"name" : "mt_chunk", "nn" : dn.attrs["nn"], "ne" : de.attrs["ne"], "eType" : eType, "nPerElem" : nPerElem, "eLabel" : eLabel})
