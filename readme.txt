Code for "Order-One Rolling Shutter Cameras": Explanations and Instructions.

The code files in this folder compute a complete classification of all minimal
problems for the relative-pose reconstruction problem with linear order-one
rolling shutter cameras. The 2D images are assumed to be configurations of
points and lines without occlusion. The camera model assumes constant rotation,
linear motion, and a rolling shutter. The order-one assumption is assumed for
cameras, and is equivalent to the camera movement being parallel to the image
plane.

The files in this folder accomplish the following workflow:

  1. Find all balanced problems

  File:   "dimension-count-infty.py"
  Input:  none.
  Output: "model-parameters-infty.m2", a list of balanced problems written as
          signature tuples.
  
  2. Compute the picture-taking map of each balanced problem and check which of
     the balanced problems are minimal. Generate, for each minimal problem,   
     Julia code in order to numerically compute the degree of the problem.

  File:   "camera-models-infty.m2"
  Input:  "model-parameters-infty.m2"
  Output: the folder "Minimal Problems" and the Julia files within. Each file
          computes the degree of a single minimal problem.

  3. Compute a lower bound for the degree of a minimal problem using monodromy 
     and interval certification.
   
  File:   every file in the Minimal Problems folder
  Input:  none
  Output: console output: calculated degree of the minimal problem.

  Note:   It is suggested to run these computations in the Julia interactive 
          interface in order to have more control over the computation. If the
          computation takes too long it can be interrupted with Ctrl-C and the
          resulting degree is still displayed and certified, which gives a lower 
          bound on the true degree.

Typical workflow from a unix terminal:

python dimension-count-infty.py
M2 camera-models-infty.m2
cd Minimal\ Problems
julia
include("J100120221.jl")


