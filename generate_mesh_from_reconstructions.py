from dolfin import *
import mshr
import sys, glob, os
import argparse
import numpy as np



#### This script was written by Benjamin Kehlet, and modified by Karl Erik Holter.

"""Script for generating DOLFIN mesh from the ECS reconstructions of Kinney, Justin P., et al. (2013).
ECS reconstructions Can be downloaded at http://www.mcell.cnl.salk.edu/models/hippocampus-ecs-2012-1/reconstructions.tar.bz2
To use this script, download the ECS reconstructions and pass the folder name of one of them as command line argument."""

dolfin.set_log_level(dolfin.TRACE)
def parse_mesh_file(filename) :
    print "Processing: ", filename
    vertices = []
    facets = []
    with open(filename, 'r') as f :
        for line in f :
            splitted = line.split()
            if splitted[0] == "Face" :
                assert len(splitted) == 5
                facets.append((int(splitted[2])-1,
                               int(splitted[3])-1,
                               int(splitted[4])-1))
            elif splitted[0] == "Vertex" :
                assert len(splitted) == 5
                vertices.append((float(splitted[2]),
                                 float(splitted[3]),
                                 float(splitted[4])))
    return mshr.CSGCGALDomain3D(np.array(vertices),
                                np.array(facets))

def get_cellular(data_path) :
    d_files = glob.glob(os.path.join(data_path, "d*.mesh"))

    d_geometry = mshr.CSGCGALDomain3D()
    for filename in d_files :
        g = parse_mesh_file(filename)
        d_geometry.insert(g)
    d_geometry.save_off("d.off")

    a_files = glob.glob(os.path.join(args.data_path, "a*.mesh"))
    a_geometry = mshr.CSGCGALDomain3D()
    for filename in a_files :
        g = parse_mesh_file(filename)
        a_geometry.insert(g)
    a_geometry.save_off("a.off")

    g_geometry = parse_mesh_file(os.path.join(args.data_path, "g.mesh"))
    g_geometry.save_off("g.off")

    cellular = d_geometry
    cellular.insert(a_geometry)
    cellular.insert(g_geometry)
    print cellular.volume()
    cellular.save_off("cellular.off")
    return cellular

def get_extra_cellular(cellular, relative_width) :
    ## full domain
    XMIN, YMIN, ZMIN = 3060, 3310, 3500
    XMAX, YMAX, ZMAX = 7130, 7000, 7580

    ## small piece of the domain
    # XMIN, YMIN, ZMIN = 3060, 3310, 3500
    # XMAX, YMAX, ZMAX = 5095, 5155, 5540 

    box = mshr.Box(Point(XMIN, YMIN, ZMIN),
                   Point(XMAX, YMAX, ZMAX)) 
    extra_cellular = mshr.CSGCGALDomain3D(box-cellular)

    print("Components: {}".format(extra_cellular.num_disconnected_components()))
    extra_cellular.keep_largest_components(1)
    extra_cellular.save_off("extra-cellular-domain.off")

    return extra_cellular

parser = argparse.ArgumentParser("parse reconstruction files, create domain and generate mesh")
parser.add_argument("data_path", metavar="data-path", type=str, action="store", help="path to folder with .mesh files")
#parser.add_argument("--check-self-intersections", default=False, action="store_true", help="")
parser.add_argument("--mesh-resolution", type=float, default=150., help="Cell resolution in mesh. Larger value gives more cells")
args = parser.parse_args(sys.argv[1:])

generator = mshr.CSGCGALMeshGenerator3D() # if this causes trouble, try mshr.TetgenMeshGenerator3D() instead

generator.parameters["detect_sharp_features"] = False
generator.parameters["mesh_resolution"] = args.mesh_resolution
m = generator.generate(get_extra_cellular(get_cellular(args.data_path), .8))
f = XDMFFile(mpi_comm_world(), "mesh.xdmf")
f.write(m)
