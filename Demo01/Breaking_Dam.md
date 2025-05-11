# OpenFOAM Breaking Dam

Reference: https://doc.cfd.direct/openfoam/user-guide-v12/dambreak#x6-230002.2

# 1 - Breaking of a Dam

Dam Break Simulation Explanation (Key Concepts)
1.	Problem Type:
    - This is a 2D simulation of a dam break.
    - It models two fluids: water and air, with a sharp interface between them (this means there’s a clear separation between water and air).
2.	Solver Used:
    - The incompressibleVoF solver is used for the simulation.
    - VoF stands for Volume of Fluid. It tracks how the two fluids (water and air) interact and evolve in each computational cell.
3.	Phase Transport Equation:
    - A phase transport equation is used to calculate the volume fraction of water and air in each cell.
    - The phase fraction is a value between 0 and 1:
    - 0 means the cell is fully air.
    - 1 means the cell is fully water.
    - For intermediate values (e.g., 0.5), the cell contains a mix of water and air.
    - This fraction tells the solver how to compute properties like density and viscosity based on the fraction of each phase in the cell.
4.	Interface Representation:
    - The interface between the two fluids (water and air) is not explicitly computed. Instead, it emerges as a result of the phase fraction.
    - This means the interface between water and air isn’t perfectly defined but is spread across a region around the actual sharp boundary (i.e., the interface isn’t a straight line but a band).
5.	Test Setup:
    - There is a column of water at rest on the left side of the tank, behind a membrane.
    - At time 0, the membrane is removed.
    - The water collapses and starts flowing down.
    - As it flows, the water impacts an obstacle at the bottom of the tank.
    - This creates complex flow patterns and air pockets (bubbles) in the water.
6.	Flow Behavior:
    - As the water collapses, the movement creates complicated flow structures and air pockets.
    - The sharp interface (where water and air meet) evolves over time, leading to a complex interaction of phases.

Visualizing the Problem:
* Initial Setup: Imagine a tank with water behind a membrane. Once the membrane is removed, the water flows down, hits the bottom, and interacts with air, creating bubbles and complex patterns.
* The solver models how water and air evolve, and how their interface (the boundary between them) changes over time.

⸻

So, the main concepts are:
* Two-phase flow (water and air).
* Volume of Fluid (VoF) method.
* Phase fraction in each computational cell.
* The interface emerges dynamically as the simulation progresses.

This setup helps simulate fluid dynamics in scenarios like dam breaks, where two fluids interact dynamically.

## 1.1 - Mesh Generation

```shell
cd $FOAM_RUN
cp -r $FOAM_TUTORIALS/incompressibleVoF/damBreakLaminar $FOAM_RUN
```

then move to ***damBreakLaminar*** director

```shell
cd damBreakLaminar
```

Once in the case directory, you can check if the blockMeshDict file exists in the system folder:

```shell
ls system/blockMeshDict
```

Now, we can generate the mesh using the blockMesh utility. Run the following command in the terminal:

```shell
blockMesh
```

Check if the blockMeshDict is correct as in the document:

```shell
nano system/blockMeshDict
```

When you run blockMesh, OpenFOAM creates the mesh in a directory called polyMesh, which is located in the constant directory of your case. This polyMesh directory contains several files that describe the mesh structure, including the mesh’s vertices, faces, and cell connections.

You can list the files in the polyMesh directory to verify

```shell
ls constant/polyMesh
```

You should see files such as ***points***, ***faces***, ***owner***, ***neighbour***, and ***boundary***. If everything is correct, the mesh has been generated successfully.


* Inside this directory, you should see several important files:
	- faces: Defines the faces of each cell.
	- points: Lists the coordinates of each vertex.
	- owner: Lists which cell is the owner for each face.
	- neighbour: Lists which neighboring cells are connected to each face.
	- boundary: Describes the boundary conditions for the simulation (i.e., the patches of the mesh like walls, inlets, outlets, etc.).


## 1.2 - Boundary Condition
The boundary ﬁle can be read and understood by the user. The user should take a look at its contents, either by opening it in a ﬁle editor or printing out in the terminal window using the cat utility

```shell
cat constant/polyMesh/boundary
```

The file might look something like this (based on the dam break tutorial you’re following):

```text
boundary
(
    leftWall
    {
        type wall;
        faces
        (
            (0 12 16 4)
            (4 16 20 8)
        );
    }

    rightWall
    {
        type wall;
        faces
        (
            (7 19 15 3)
            (11 23 19 7)
        );
    }

    lowerWall
    {
        type wall;
        faces
        (
            (0 1 13 12)
            (1 5 17 13)
            (5 6 18 17)
            (2 14 18 6)
            (2 3 15 14)
        );
    }

    atmosphere
    {
        type patch;
        faces
        (
            (8 20 21 9)
            (9 21 22 10)
            (10 22 23 11)
        );
    }
);
```

1. Five boundary patches are defined:
- leftWall
- rightWall
- lowerWall
- atmosphere
- defaultFaces

2. Atmosphere Patch
- This is a standard patch with no special attributes.
- It simply allows you to apply boundary conditions for outflow or inflow of fluid.

3. defaultFaces Patch
- This patch represents faces omitted from the mesh (from the blockMeshDict).
- It’s typically marked as empty for 2D problems, meaning no flow is solved in the third direction.

4. Wall Patches (leftWall, rightWall, lowerWall)
- These are marked as wall type patches.
- A wall does not have geometric or topological information but serves to identify the boundary.
- It’s important for certain models, such as turbulence models, to calculate things like the distance to the nearest wall.

5. Surface Tension and Wall Adhesion (VoF Method)
- VoF (Volume of Fluid) uses a sharp interface between fluids.
- Surface tension models can include wall adhesion at the interface between fluid and wall.
- alphaContactAngle boundary condition can be used for this.
- This example ignores surface tension, using a simpler approach by applying zeroGradient condition to the alpha (phase fraction) field on the walls.

6. Boundary Conditions on the Top Boundary (Atmosphere)
- The top boundary is open to the atmosphere, so it needs to allow both inflow and outflow.
- Two types of boundary conditions are applied:
- prghTotalPressure: Applied to the pressure field minus the hydrostatic component.
- pressureInletOutletVelocity: Applied to the velocity field, sets zeroGradient unless there is inflow, in which case it uses fixedValue.

7. Pressure Boundary Conditions at Walls
- On all wall boundaries, the fixedFluxPressure boundary condition is used.
- This adjusts the pressure gradient to match the velocity boundary condition, especially for gravity and surface tension effects.

8. defaultFaces Patch for 2D
- The defaultFaces patch, which represents the front and back planes in 2D simulations, is usually set to empty.

## 1.3 - Phases

the fluidPhase is specified in the ***phaseProperties*** file within the ***constant*** directory. To see this information:

```shell
cat constant/phaseProperties
```

```text
16

17  phases          (water air);

18

19  sigma           0.07;

20

21

22// ************************************************************************* //
```

You should see something like this.

The file defines two phases: water and air. The phase fraction (how much of each phase is present) is calculated for the listed phases, but not for the last phase (in this case, air). Since there are only two phases, only the equation for the water phase fraction (called alpha.water) is solved, and it is found in the 0 directory.

Additionally, the phaseProperties file includes information about the surface tension between the two phases, represented by sigma, measured in N/m.

## 1.4 - Setting Initial Fields

In this step, we set the initial values for the phase fraction, which determines how the water and air phases are distributed across the computational domain.

### 1. Phase Fraction
The **phase fraction** represents how much of each phase (water and air) is present in each cell. The value of the phase fraction can range from 0 to 1:
- **1** means the entire cell is filled with water (water phase).
- **0** means the entire cell is filled with air (air phase).

### 2. Default Field Values
We start by defining **default** values for the fields. These fields represent variables we will work with in the simulation, such as the phase fraction.
- For example, `alpha.water` is set to **0** initially. This means the domain is assumed to have no water and only air at the start.

### 3. Execute
```shell
setFields
```

If you wish to see what is inside the setField, you can see by:

```shell
nano system/setFieldsDict
```

then you should have something like this:

```
/*--------------------------------*- C++ -*----------------------------------*\
 |                                                                           |
 | OpenFOAM: The Open Source CFD Toolbox                                      |
 | Website: https://openfoam.org                                               |
 | Version: 12                                                                |
 |---------------------------------------------------------------------------|
\*---------------------------------------------------------------------------*/

defaultFieldValues
(
    volScalarFieldValue alpha.water 0
    volScalarFieldValue tracer.water 0
    volScalarFieldValue tracer.air 0
);

regions
(
    boxToCell
    {
        box (0 0 -1) (0.1461 0.292 1);
        fieldValues
        (
            volScalarFieldValue alpha.water 1
        );
    }

    boxToCell
    {
        box (0 0 -1) (0.1461 0.146 1);
        fieldValues
        (
            volScalarFieldValue tracer.water 1
        );
    }

    boxToCell
    {
        box (0 0.292 -1) (0.1461 0.438 1);
        fieldValues
        (
            volScalarFieldValue tracer.air 1
        );
    }
);
```

## 1.5 - Fluid Properties

In OpenFOAM, we define the physical properties of the fluids (air and water) in separate files:

- `physicalProperties.air` for air
- `physicalProperties.water` for water

### Key Properties:
1. **Viscosity Model**:
   - Set to `constant`, meaning the viscosity does not change with the flow.
   
2. **Kinematic Viscosity (`nu`)**:
   - For **air**: `1.48e-05` m²/s
   - For **water**: `1e-06` m²/s

3. **Density (`rho`)**:
   - For **air**: `1 kg/m³`
   - For **water**: Typically around `1000 kg/m³` (but you can change this based on your simulation).

### Example of `physicalProperties.air` file:

```plaintext
viscosityModel  constant;
nu              1.48e-05;
rho             1;
```

-	viscosityModel constant: The viscosity of air is constant.
-	nu: Kinematic viscosity for air.
-	rho: Density of air.

If you wish to see the specificatio of ***air*** or ***water***. It is located in the ***constant*** directory. For example:

```shell
nano constant/physicalProperties.air
```

## 1.6 - Gravity

In OpenFOAM, gravitational acceleration is specified using a `uniformDimensionalVectorField`. This means the field for gravity is uniform across the entire domain.

- The gravitational acceleration (`g`) is defined in a file named `g` within the `constant` directory.
- This file contains two components:
  - **dimensions**: Defines the physical dimensions of the gravitational vector (e.g., meters per second squared).
  - **value**: Specifies the gravity vector, in this case, representing gravitational acceleration along the negative Y-axis, i.e., `(0 -9.81 0)` (m/s²).

Here’s an example of the `g` file contents:

```
dimensions [0 1 -2 0 0 0 0];   // Dimensions for acceleration in m/s²
value (0 -9.81 0);              // Gravitational value in the negative Y direction
```

## 1.7 - Turbulence Modelling

- In OpenFOAM, the turbulence modelling method can be selected during runtime using the `simulationType` keyword in the `momentumTransport` dictionary.
- In this example, turbulence is not needed, so the `simulationType` is set to `laminar`, which means no turbulence model will be applied.

Example:

```plaintext
simulationType  laminar;   // No turbulence model
```

## 1.8 - Time Step Control

The simulation of the **damBreakLaminar** case is fully transient, so the time step requires careful attention. The **Courant Number (Co)** is an important consideration when determining the time step. It is a dimensionless parameter that can be defined for each cell as:

```plaintext
Co = (δt |U|) / δx
```

Where:
- δt = Time step
- |U| = Magnitude of the velocity through that cell
- δx = Cell size in the direction of the velocity

This equation represents the relationship between the time step, velocity, and cell size.

### Key Points:
1.	Explicit Solutions:
-	For explicit solutions, stability requires the Co to be less than 1. This means the time step must be sufficiently small. Stricter limits exist depending on the choice of advection scheme.
-	If Co exceeds 1, it may result in instability or incorrect solutions, especially with advection schemes.
2.	Implicit Solutions:
-	In implicit solutions, the stability does not depend on the maximum Co, but increasing Co beyond 1 may affect temporal accuracy.
-	Implicit methods allow larger time steps, but care must be taken when Co increases significantly beyond 1 as it may affect the solution’s accuracy.

### Time Step Control with MULES:
-	The time step is particularly important in interface-capturing, where the simulation involves multiple phases with a sharp interface (such as air and water).
-	The incompressibleVoF solver uses the multidimensional universal limiter for explicit solution (MULES) to maintain boundedness of the phase fraction.
-	The original explicit version of MULES requires an upper limit of Co ≈ 0.25 for stability.
-	However, there is also a semi-implicit version of MULES specified by the MULESCorr switch in the fvSolution file. This allows for larger values of Co but is determined by the level of temporal accuracy.

### Why Adjust Time Step?

- The adjustment of the time step ensures that the simulation remains stable and accurate as the flow evolves, especially in situations involving phase change or complex interfaces (e.g., dam break simulations).

### Explanation Summary:
- The **Co** equation defines the relationship between the time step, velocity, and cell size.
- **Explicit solutions** require **Co < 1** for stability.
- **Implicit solutions** allow larger time steps but require more attention to temporal accuracy.
- **MULES** (used for phase fraction in multi-phase simulations) limits **Co** for stability.
- Use **automatic time step control** to dynamically adjust time steps during the simulation and ensure stability.
- The **controlDict** file allows you to control the writing intervals and adjust the time step automatically.


## 1.9 - Discretisation Schemes

The **MULES** method, used by the **incompressibleVoF** solver, maintains the boundedness of the phase fraction independently of the underlying numerical scheme, mesh structure, etc. The choice of convection schemes is not restricted to those that are strongly stable or bounded, such as upwind differencing.

The convection schemes settings are specified in the **divSchemes** sub-dictionary of the **fvSchemes** dictionary. Here is an example for the convection term in the momentum equation, **∇•(ρUU)**, which is denoted by the **div(rhoPhi,U)** keyword. This uses **Gauss linearUpwind grad(U)** to produce good accuracy.

Additionally, for the **∇•(Uα)** term, the **div(phi,alpha)** keyword is used with a bespoke **interfaceCompression** scheme. This scheme helps control the compression of the interface where:
- `0` corresponds to no compression.
- `1` corresponds to conservative compression.
- Anything larger than `1` enhances compression of the interface.

In this example, we use **1.0** for the coefficient, which corresponds to **conservative compression**.

### Convection Schemes Configuration Example:

```plaintext
divSchemes
{
    div(phi,alpha)  Gauss interfaceCompression vanLeer 1;
    div(rhoPhi,U)   Gauss linearUpwind grad(U);
    div(phi,k)      Gauss upwind;
    div(phi,epsilon) Gauss upwind;
    div(phi,omega)  Gauss upwind;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
```

### Explanation:
-	div(phi,alpha): This term controls the compression of the interface. We use the interfaceCompression scheme with a coefficient value of 1.0 (conservative compression).
-	div(rhoPhi,U): Convection term in the momentum equation, where Gauss linearUpwind grad(U) is used for good accuracy in the flow velocity field.
-	Other terms: These terms include div(phi,k), div(phi,epsilon), etc., and use commonly employed schemes for their discretization.

### Explanation in Simple Terms:
- **MULES** is a method that ensures the phase fraction (in this case, water and air) is calculated correctly, and it doesn't depend on the numerical scheme, mesh, etc.
- The convection schemes handle how quantities like velocity and pressure are calculated over time and space.
- We use **interfaceCompression** for handling the interface between two phases (water and air) and set it to **1.0** for conservative compression.
- You need to update the **fvSchemes** file to specify how to handle these terms and ensure good accuracy in your simulation.

### Action Required:
You only need to **check or modify the `fvSchemes` file** to ensure that these schemes are defined correctly. No additional action is needed beyond configuring the file.

## 1.10 - Linear-Solver Control

In the **fvSolution** file, the sub-dictionary in **solvers** for **alpha.water** contains elements that are specific to the **MULES** algorithm. Here is an example configuration for the solver:

```plaintext
"alpha.water.*"
{
    nAlphaCorr      2;
    nAlphaSubCycles 1;
    MULESCorr       yes;
    nLimiterIter    5;
    solver          smoothSolver;
    smoother        symGaussSeidel;
    tolerance       1e-8;
    relTol          0;
}
```

### Explanation:
-	nAlphaCorr: Specifies the number of iterations for the phase fraction equation. It is used to correct nonlinearities in the advection of the phases. A typical value is 2.
-	nAlphaSubCycles: Specifies the number of sub-cycles for solving alpha.water in each step. Typically, this is set to 1.
-	MULESCorr: If set to yes, this activates the semi-implicit MULES method, which first calculates an implicit upwind solution and then applies MULES as a higher-order correction.
-	nLimiterIter: Refers to the number of iterations for the limiter. The limiter keeps the phase fraction within the range [0, 1]. A rule of thumb is that it should be set to 3 + Co, where Co is the Courant number, to maintain boundedness.
-	solver: Defines the solver to use. In this case, we use smoothSolver which is used for smooth convergence.
-	smoother: Specifies the type of smoother used in the solver. symGaussSeidel is commonly used for smoother iterations.
-	tolerance: Defines the convergence tolerance for the solver, which is set to 1e-8 for this case.
-	relTol: This is the relative tolerance, set to 0 in this example, which means the absolute tolerance is used.

### Explanation in Simple Terms:

- This section configures the **solver** settings used to calculate the phase fraction (`alpha.water`), which is crucial in **MULES** method for handling two-phase flow (water and air).
- **MULES** is an algorithm that ensures the phase fraction stays within the 0-1 bounds.
- These settings control how many times the solver runs to correct any errors in the phase fraction calculation.
- You typically do not need to change these settings unless you are running a very specific or complex case.

## 1.11 - Run the simulation
to Run
```shell
foamRun | tee log
```

to observe
```shell
paraFoam
```
then in ***Properties*** pannel Apply, in **Coloring** section choose ***alpha.water*** to see water. In the **toolbar** click run to see the result.

