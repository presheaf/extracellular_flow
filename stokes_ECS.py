import sys, os
import numpy as np
from dolfin import *

"""Solve a Stokes problem on mesh of reconstructed ECS."""

### helper functions for running in parallel
def is_primary_process():
    """True if MPI rank is 0, otherwise False"""
    return MPI().rank(mpi_comm_world()) == 0

def print_message(message):
    """Helper function for handling output in parallel."""
    if is_primary_process():
        print message
    else:
        pass

if not len(sys.argv) == 5:
    print_message(
        "Usage: python stokes_ECS.py mesh_fn output_dir P1P1 x/y/z"
    )
    sys.exit(1)

output_dir = sys.argv[2]
if is_primary_process():
    try:
        os.makedirs(output_dir)
    except OSError, e:
        pass
    

### read input
mesh_fn = sys.argv[1]
use_P1P1 = (sys.argv[3] == "P1P1")  # P1P1 for stabilized P1P1, otherwise P2P1
flow_direction = sys.argv[4].lower()

if flow_direction not in ["x", "y", "z"]:
    print_message("Please use x, y or z for flow direction.")
    sys.exit(1)

print_message("Using flow direction " + flow_direction)


class ReconstructedECSDomain(object):
    """Class for loading ECS mesh from file and marking boundary."""
    def __init__(self, fn):
        self.filename = fn
    
    def partition_boundary(self):
        """Returns a FacetFunction labeling the subdomains with ints."""
        boundary_parts = FacetFunction("size_t", self.mesh)
        
        xmin, xmax, ymin, ymax, zmin, zmax = self.get_coord_minmax(self.mesh)

        near_eps = 1E-3
        
        class HoleBoundary(SubDomain): 
            def inside(self, x, on_boundary):
                # hard to specify hole boundaries, so mark everything
                # and then re-mark everything which is not a hole
                return on_boundary
        HoleBoundary().mark(boundary_parts, 1)

        class LeftBoundary(SubDomain): 
            def inside(self, x, on_boundary):
                return (near(x[0], xmin, eps=near_eps) and on_boundary)

        LeftBoundary().mark(boundary_parts, 2)

        class RightBoundary(SubDomain): 
            def inside(self, x, on_boundary):
                return (near(x[0], xmax, eps=near_eps) and on_boundary)

        RightBoundary().mark(boundary_parts, 3)

        
        class BackBoundary(SubDomain):
            def inside(self, x, on_boundary):
                return (near(x[1], ymax, eps=near_eps) and on_boundary)

        BackBoundary().mark(boundary_parts, 4)
        
        class FrontBoundary(SubDomain):
            def inside(self, x, on_boundary):
                return (near(x[1], ymin, eps=near_eps) and on_boundary)

        FrontBoundary().mark(boundary_parts, 5)


        class TopBoundary(SubDomain):
            def inside(self, x, on_boundary):
                return (near(x[2], zmax, eps=near_eps) and on_boundary)

        TopBoundary().mark(boundary_parts, 6)


        class BotBoundary(SubDomain):
            def inside(self, x, on_boundary):
                return (near(x[2], zmin, eps=near_eps) and on_boundary)

        BotBoundary().mark(boundary_parts, 7)

        self._subdomain_data = boundary_parts
        return self._subdomain_data


    def get_coord_minmax(self, mesh):
        """Return list of xmin, xmax, ymin, ymax, zmin, zmax of a mesh."""
        coords = mesh.coordinates()
        xcoords, ycoords, zcoords= coords[:, 0], coords[:, 1], coords[:, 2]
        xmin = min(xcoords)
        xmax = max(xcoords)
        ymin = min(ycoords)
        ymax = max(ycoords)
        zmin = min(zcoords)
        zmax = max(zcoords)

        return [xmin, xmax, ymin, ymax, zmin, zmax]
        

    def load_mesh(self):
        """Load mesh from file."""
        if self.filename[-5:] == ".xdmf":
            mesh = Mesh()
            f = XDMFFile(mpi_comm_world(), self.filename)
            f.read(mesh, True) # True = try to load mesh partition from file
            return mesh
        else:
            return Mesh(self.filename) 

    @property
    def mesh(self):
        """Mesh of domain as DOLFIN object."""
        try:
            return self._mesh
        except AttributeError:
            print_message("Starting mesh creation.")
            self._mesh = self.load_mesh()
            return self._mesh

    @property
    def subdomain_data(self):
        """FacetFunction labeling the subdomains with ints."""
        try:
            return self._subdomain_data
        
        except AttributeError:
            print_message("Paritioning boundary.")
            self._subdomain_data = self.partition_boundary()
            return self._subdomain_data
    

    def measure(self, measure_type):
        """Return DOLFIN ds/dx objects."""
        return Measure(measure_type, domain=self.mesh,
                       subdomain_data=self.subdomain_data)

    def refine(self):
        """Subdivide every tet of mesh."""
        refined_mesh = refine(self.mesh)
        del self._mesh
        if hasattr(self, "_subdomain_data"):
            del self._subdomain_data
        self._mesh = refined_mesh

    
class StokesProblem(object):
    """Class for solving a Stokes problem on a given domain."""
    def __init__(self, argdict, domain,
                 noslip_subdomains, p_N_subdomains, symm_subdomains):
        self.mu = argdict["mu"]

        self.atol = argdict["atol"]
        self.rtol = argdict["rtol"]
        self.gmres_restarts = 50
        self.max_its = 30000
        
        self.initial_guess = argdict["initial_guess"]
        self.print_convergence_output = argdict["print_convergence_output"]

        self.solver = argdict["solver"]
        self.preconditioner = argdict["preconditioner"]
        
        self.domain = domain
        self.beta = argdict["beta"]

        self.use_P1P1 = argdict["use_P1P1"]
        self.W = self.function_space()
        self._matrices = None

        self.u_noslip_subdomains = noslip_subdomains
        self.p_neumann_subdomains = p_N_subdomains
        self.symmetry_subdomains = symm_subdomains

    @staticmethod        
    def default_arg_dict():
        """Dict of necessary parameters as well as default values."""
        return  {
            "mu": 0.8E-3,                     # Viscosity of fluid.
            "atol": 1E-5, #should be lowered
            "rtol": 1E-5,
            "initial_guess": None,
            "print_convergence_output": True,
            "preconditioner": "hypre_amg",
            "solver": "gmres",
            "beta": 0.2,                    # Stabilization parameter.
            "use_P1P1": True,
            
        }


    def function_space(self):
        """Create a P1P1 or a P2P1 mixed FunctionSpace."""
        elt = self.domain.mesh.ufl_cell()

        if self.use_P1P1:
            print_message("Creating P1-P1 FunctionSpace.")
            P1_u = VectorElement("Lagrange", elt, 1)
            P1_p = FiniteElement("Lagrange", elt, 1)
            P1P1 = MixedElement([P1_u, P1_p])
            return FunctionSpace(self.domain.mesh, P1P1)
        else:
            print_message("Creating P2-P1 FunctionSpace.")
            P2 = VectorElement("Lagrange", elt, 2)
            P1 = FiniteElement("Lagrange", elt, 1)
            P2P1 = MixedElement([P2, P1])
            return FunctionSpace(self.domain.mesh, P2P1)

        
    def u_noslip_bcs(self, subdomains=[]):
        """Supply list of subdomains, return list of noslip DirichletBCs"""
        return [DirichletBC(self.W.sub(0), Constant((0, 0, 0)),
                            self.domain.subdomain_data, subdomain_id)
                for subdomain_id in subdomains]

    
    def symmetry_bcs(self, subdomains={}):
        """Returns a list of Dirichlet symmetry bcs on the specified 
        subdomains. subdomains should be a dict {i: dim} where i is the
        subdomain id and dim is the coordinate vector perpendicular to subdomain i
        (0 if subdomain i is perpendicular to the x-axis, 1 if y, 2 if z)"""

        return [DirichletBC(self.W.sub(0).sub(subdomains[subdomain_id]),
                            Constant(0), self.domain.subdomain_data,
                            subdomain_id)
                for subdomain_id in subdomains]

    
    def solve(self):
        ## name a bunch of things
        ds = self.domain.measure("ds")
        u, p = TrialFunctions(self.W)
        v, q = TestFunctions(self.W)
        w = Function(self.W)
        n = FacetNormal(self.domain.mesh)
        beta = self.beta
        h = CellSize(self.domain.mesh)

        ## define variational problem
        mu = self.mu
        if self.use_P1P1:
            print_message("Using stabilized variational form with beta = {}".format(beta))
            a = (mu*inner(grad(u), grad(v))*dx + div(v)*p*dx + div(u)*q*dx
                 - beta*h*h*inner(grad(p), grad(q))*dx)
            b = (mu*inner(grad(u), grad(v))*dx + p*q/mu*dx
                 + beta*h*h*inner(grad(p), grad(q))*dx)
        else:
            print_message("Using unstabilized variational form")
            a = mu*inner(grad(u), grad(v))*dx + div(v)*p*dx + div(u)*q*dx
            b = mu*inner(grad(u), grad(v))*dx + p*q/mu*dx

        ## create BCs
        print_message("Creating noslip Dirichlet BCs for velocity.")
        dirichlet_bcs = self.u_noslip_bcs(self.u_noslip_subdomains)

        if len(self.symmetry_subdomains) > 0:
            print_message("Creating symmetric Dirichlet BCs for velocity.")
            symmetry_bcs = self.symmetry_bcs(self.symmetry_subdomains)
            dirichlet_bcs.extend(symmetry_bcs)
        
        print_message("Creating Neumann BCs.")
        L = sum([
            neumann_bc * dot(v, n) * ds(subdomain)
            for subdomain, neumann_bc in self.p_neumann_subdomains.iteritems()
        ])

        
        ## build linear system
        print_message("Assembling matrices.")
        (A, bb) = assemble_system(a, L, dirichlet_bcs)
        (P, _) = assemble_system(b, L, dirichlet_bcs)
        self._matrices = (A, P)


        ## set solver parameters
        solver = PETScKrylovSolver(self.solver, self.preconditioner)
        if self.solver == "gmres":
            solver.ksp().setGMRESRestart(self.gmres_restarts)
        solver.set_operators(A, P)
        
        kpars = parameters["krylov_solver"]
        kpars["absolute_tolerance"] = self.atol # default: 1E-15
        kpars["relative_tolerance"] = self.rtol # default: 1E-9
        kpars["maximum_iterations"] = self.max_its # default: ?

        if self.print_convergence_output:
            kpars["report"] = True
            kpars["monitor_convergence"] = True
            
        if self.initial_guess is not None:
            print_message("Using nonzero initial guess.")
            kpars["nonzero_initial_guess"] = True
            w2 = interpolate(self.initial_guess, self.W)
            w.vector()[:] = w2.vector()[:]
            
        solver.update_parameters(kpars)

        
        ## solve linear system
        print_message("Starting solve.")
        solver.ksp().setConvergenceHistory()
        num_its = solver.solve(w.vector(), bb)
        self.residuals = solver.ksp().getConvergenceHistory()
        print_message("Solving complete in {} iterations.".format(num_its))

        u, p = w.split()
        
        return u, p


## handle parameters 
solver = "gmres"
preconditioner = "hypre_amg"
print_message("Using solver {} and preconditioner {}".format(solver, preconditioner))

p_argdict = StokesProblem.default_arg_dict()
p_argdict.update(
    {
        "atol": 1E-7, "rtol": 1E-7,
        "print_convergence_output": True,
        "mu": 0.8E-3,           
        "solver": solver,
        "preconditioner": preconditioner,
        "beta": 0.2,
        "use_P1P1": use_P1P1
    }
)

ppx = 133E-6 # pressure differential of 1 mmHg/mm (in Pa/nanometer)

## create symmetric and no-slip BCs in proper flow direction.
## for performance, use a linear function as initial guess for pressure
symm_bcs = {
    2: 0, 3: 0,
    4: 1, 5: 1,
    6: 2, 7: 2
}
dirichlet_bcs = [1]

if flow_direction == "x":
    del symm_bcs[2], symm_bcs[3]
    initial_guess = Expression(("0", "0", "0",
                                "x[0] * dp"), dp=ppx, degree=1)
    p_bc = initial_guess[3]
    neumann_bcs = {i: p_bc for i in [2, 3]} 

    
elif flow_direction == "y":
    del symm_bcs[4], symm_bcs[5]
    initial_guess = Expression(("0", "0", "0",
                                "x[1] * dp"), dp=ppx, degree=1)
    p_bc = initial_guess[3]
    neumann_bcs = {i: p_bc for i in [4, 5]} 
    
else:                           # assume == "z"
    del symm_bcs[6], symm_bcs[7]
    initial_guess = Expression(("0", "0", "0",
                                "x[2] * dp"), dp=ppx, degree=1)
    p_bc = initial_guess[3]
    neumann_bcs = {i: p_bc for i in [6, 7]} 



p_argdict["initial_guess"] = initial_guess


### solve problem
domain = ReconstructedECSDomain(mesh_fn)
problem = StokesProblem(p_argdict, domain,
                        dirichlet_bcs, neumann_bcs, symm_bcs)

u, p = problem.solve()


### save output to file
if is_primary_process():
    np.save("{}/convergence_residuals".format(output_dir), problem.residuals)

f = File("{}/u.pvd".format(output_dir))
f << u

g = File("{}/p.pvd".format(output_dir))
g << p
