M2
restart

-- computations with the finite list of candidate minimal problems
-- for rolling shutter cameras with linear movement and no rotation.

-- mod out the action of the Euclidean group (6) + global scaling (1)
-- by fixing 7 of 8 camera parameters of the first camera.

-- dramatis personae:
  -- Z:      camera center as a point in space
  -- Rot:    (fixed) rotation of the camera using cayley coordinates
  -- DTilde: movement vector of the camera after standardization
  -- A:      free points
  -- B:      third collinear points
  -- C:      fourth collinear points
  -- D:      free lines
  -- E:      lines incident to one point
  -- F:      lines incident to two points
  -- pInfty: point at infinity associated to a particular camera 
  -- incid:  number of lines incident to A_0
  -- extra:  1 or 0 depending on whether we have given pInfty or not.
  -- imgA:   image point of A in 2d
  -- imgB:   similar
  -- imgC:   similar
  -- imgD:   image conic of D in 2d
  -- imgE:   similar
  -- imgF:   similar

-- create variable 'parameters' that encodes the minimal problems.
load "model-parameters-infty.m2"

-- an iterator. iter(v) gives the next v. set iterNext in the code.
iter = v -> (
  iterNext = iterNext+1;
  v_(iterNext-1)
)

-- normalized cayley rotation given three parameters.
cayleyRotation = v -> (
  w := 1; -- coordinate at infinity
  x := v_0;
  y := v_1;
  z := v_2;
  res1 := {w^2 + x^2 - y^2 - z^2, 2*x*y - 2*z*w, 2*y*w + 2*x*z};
  res2 := {2*z*w + 2*x*y, w^2 - x^2 + y^2 - z^2, 2*y*z - 2*x*w};
  res3 := {2*x*z - 2*y*w, 2*x*w + 2*y*z, w^2 - x^2 - y^2 + z^2};
  matrix{res1, res2, res3} -- corresponds to quaternion w + ix + jy +kz
)

-- take a picture of point X from the standard camera centered at Z.
-- this function is documentation and not supposed to run.

-- takePictureNoRot = (Z, X) -> (
--  xTilde := (X-Z);
--  vector{xTilde_0/xTilde_2, xTilde_1/xTilde_2}
-- )
 
-- infer time t from X and movement vector D of standard camera [id|0].
-- formula solves "takePictureNoRot(X, t*D)_0 = t" for t.
timeNoRot = (X, D) -> (
  (X_0)/(X_2 + D_0)
)

-- use the inferred time t to take a picture of the point X
-- with camera parameters Z, Rot, DTilde, r0.
-- DTilde is the direction of movement of the camera center, *after*
-- having transformed the camera into a standard camera [id|0]
takePic = (Cam, X) -> (
  (Z, Rot, DTilde) := Cam;
  X = sub(X, frac ring X);
  DTilde = sub(DTilde, frac ring DTilde);
  Rot = sub(Rot, frac ring Rot);
  xTransf := cayleyRotation(Rot) * (X - Z);
  timeInferred := timeNoRot(xTransf, DTilde);
  xTilde := (xTransf - timeInferred * DTilde);
  vector{xTilde_0/xTilde_2, xTilde_1/xTilde_2}
)

-- calculate the image conic of a line passing through P and Q.
-- the conic equation is -v x^2 + xy + c_02 x + c_12 y + c_22.
-- this is normalized with c_01 = 1 and pInfty = (1, v).

takePicLine = (Cam, P, Q) -> (
  origRing = ring P;
  T = QQ[gens origRing, t, s_0, s_1];
  PT = sub(P, T);
  QT = sub(Q, T);
  CamT = toSequence(apply(toList(Cam), u -> sub(u, T)));
  L = (1-t)*PT + t*QT;
  tP = takePic(CamT, L);
  myMap = {numerator(s_0 - tP_0), numerator(s_1 - tP_1)};
  myIdeal = eliminate(t, ideal(myMap));
  myGen = sub(myIdeal_0, QQ[gens origRing, t][s_0, s_1]);
  resCoeffs = (coefficients (myGen, Monomials => {s_0, s_1, 1}))_1_0;
  resCoeffsOrig = apply(resCoeffs, u -> sub(u, origRing));
  c01 = (coefficients (myGen, Monomials => {s_0*s_1}))_1_0_0;
  c01Orig = sub(c01, origRing);
  for i in 0..2 list resCoeffsOrig_i/c01Orig
)

-- calculate the coordinates of the rational map phi associated to
-- the input parameter.
mapCoords = parameter -> (
  
  (a, b, c, d, e, f, m, incid, extra) = parameter;
  
  -- create source and target rings for the taking pictures map
  rDim = 3*a + b + c + 4*d + 2*e + 8*m - 7;
  dimensionCheck = (rDim == m*(extra + 2*a + 2*b + c + 3*d + 2*e + f));
  R = QQ[x_1..x_rDim];
  
  -- number of coordinates of the map. may exceed rDim if C points are
  -- present.
  excessTargetDim = m*(extra + 2*a + 2*b + 2*c + 3*d + 2*e + 1*f);
  S = QQ[y_1..y_excessTargetDim]; -- target ring not being used yet.
      
  iterNext = 1; -- keep track of current variable. now it's x_1.
  
  -- parameters of the first camera
  Z_0 = sub(vector{0, 0, 0}, R);
  Rot_0 = sub(vector{0, 0, 0}, R);
  DTilde_0 = vector{1, iter(x), 0};
  Cam_0 = (Z_0, Rot_0, DTilde_0);
  
  --parameters of the other cameras
  for j in 1..(m-1) do (
    Z_j = vector{iter(x), iter(x), iter(x)};
    Rot_j = vector{iter(x), iter(x), iter(x)};
    DTilde_j = vector{iter(x), iter(x), 0};
    Cam_j = (Z_j, Rot_j, DTilde_j);
  );
  
  -- A-F: point-line configuration in space.
  for i in 0..(a-1) do (
    A_i = vector{iter(x), iter(x), iter(x)};
  );
  
  -- have b <= 2. attach B_0 to A_0 and A_1; attach B_1 to A_0 and A_2
  -- this suffices since if b = 2 then a = 3, thus two of the B
  -- must share an A.
  for i in 0..(b-1) do (
    B_i = A_0 + iter(x) * (A_0 - A_(i+1));
  );
  
  -- if c > 0 then b = 1. Thus all C attach to A_0 and A_1.
  for i in 0..(c-1) do (
    C_i = A_0 + iter(x) * (A_0 - A_1);
  );
  
  -- represent a line as a pair of affine points it passes through. 
  for i in 0..(d-1) do (
    D_i = {vector{1, iter(x), iter(x)},
           vector{0, iter(x), iter(x)}};
  );
  
  -- we only have three combinatorially ambiguous balanced problems:
  -- (2, 0, 0, 0, 2, 0), (3, 0, 0, 0, 1, 1) and (3, 0, 0, 0, 2, 0, 2)
  -- the parameter incid, the number of lines passing through A_0,
  -- resolves the ambiguity. Now build the E and F lines with this
  -- in mind. To do this: create variables
  -- alpha, the number of E lines passing through A_0, and beta, the
  -- number of F lines passing through A_0. Can do this by inspection
  -- of the balanced problems. There's only one case with both
  -- 'E' and 'F' lines. That case we can divide into beta=1 and beta=0
  
  alpha = min(e, incid); -- will equal incid unless e=1,f=1,incid=2.
  beta = incid - e; -- will equal incid unless e=1,f=1,incid=1 or 2.
  
  -- we have e <= 3. thus alpha completely determines the
  -- configuration of the E. Wlog E can always attach to 'A' points
  -- since if b > 0 then e <= a.
  for i in 0..(alpha-1) do (
    E_i = {A_0, vector{0, iter(x), iter(x)}};
  );

  for i in alpha..(e-1) do (
    E_i = {A_(i-alpha+1), vector{0, iter(x), iter(x)}};
  );

  -- we have f <= 3. the +b ensures we are not reusing 'B' lines
  for i in 0..(beta-1) do (
    F_i = {A_0, A_(i+1+b)}; 
  );

  -- had to reverse the order to make 5 point cases correct in earlier
  -- version.
  for i in beta..(f-1) do (
    F_i = {A_(a-i+beta-2), A_(a-1)};
  );
  
  coords = {};
  
  for j in 0..(m-1) do (
    
    --only record pInfty if it is given
    if (extra == 1) then
	    imgPInfty_j = (DTilde_j)_1/(DTilde_j)_0;
    
    imgA_j = for i in 0..(a-1) list takePic(Cam_j, A_i);
    
    imgB_j = for i in 0..(b-1) list takePic(Cam_j, B_i);
    
    -- images of C-points are identified by their projection to the
    -- x-axis.
    imgC_j = for i in 0..(c-1) list ((takePic(Cam_j, C_i))_0);
    
    -- parametrize D conics by c02, c12, c22
    coeffsD_j = for i in 0..(d-1) list takePicLine(Cam_j, D_i_0, D_i_1);
    imgD_j = for i in 0..(d-1) list {coeffsD_j_i_0, coeffsD_j_i_1,
                                     coeffsD_j_i_2};
    
    -- parametrize E conics by c02, c12
    coeffsE_j = for i in 0..(e-1) list takePicLine(Cam_j, E_i_0, E_i_1);
    imgE_j = for i in 0..(e-1) list {coeffsE_j_i_0, coeffsE_j_i_1};  
    
    -- parametrize F conics by c02
    coeffsF_j = for i in 0..(f-1) list takePicLine(Cam_j, F_i_0, F_i_1);
    imgF_j = for i in 0..(f-1) list {coeffsF_j_i_0};

    -- build coordinates of the phi map.
    -- only record pInfty if it is given
    if (extra == 1) then
        coords = coords | {imgPInfty_j};
    coords = coords | (flatten apply(imgA_j, u -> entries u));
    coords = coords | (flatten apply(imgB_j, u -> entries u));
    coords = coords | (flatten apply(imgC_j, u -> u));
    coords = coords | (flatten imgD_j);
    coords = coords | (flatten imgE_j);
    coords = coords | (flatten imgF_j);
  );
  
  return coords;
)

-- evaluate a rat. function at the point p, except for argument i.
evalPartial = (p,i,f) -> (
  args := gens ring numerator f;
  pHat := insert(i, args_i, drop(p,{i,i})); -- replace p_i by x_i
  num := sub(numerator(f), matrix{pHat});
  denom := sub(denominator(f), matrix{pHat});
  num/denom
)

-- compute the derivative of a rat. function f by the variable x.
fracDiff = (x,f) -> (
  g := numerator f;
  h := denominator f;
  g' := diff(x,g);
  h' := diff(x,h);
  (g'*h - g*h')/h^2
)

-- evaluate a *univariate* rat. function at a point p.
evalFrac = (p,f) -> (
  num := sub(numerator(f), matrix{p});
  denom := sub(denominator(f), matrix{p});
  num/denom
) 

-- computes the Jacobian of a rational function f at the point p.
-- assumes f is given as a list of rational functions in the same ring.
evalJac = (p,f) -> (
  args = gens ring numerator f_0;
  for j in 0..(length(f)-1) do (
    for i in 0..(length(args)-1) do (
    partial := evalPartial(p,i,f_j);
      partialDiff := fracDiff(args_i, partial);
      jac_(j,i) = evalFrac(p,partialDiff);
    );
  );
  matrix(
    for j in 0..(length(f)-1) list (
      for i in 0..(length(args)-1) list jac_(j,i)
    )
  )
)

-- check balanced problems for minimality by computing the jacobian
-- at a random point. This loop is better executed in the interactive
-- interface, since the random selection of a point might cause a
-- division by zero. In that case, simply re-run the loop. For this
-- reason, the loop is commented out, but can be commented in again
-- for testing in the interactive interface.

--for param in parameters do (
--  phiCoords = mapCoords param;
--  -- try to avoid dividing by zero
--  try (
--  p = for i in 0..(rDim-1) list (random QQ);
--  myJac = evalJac(p,phiCoords);
--  ) then () else (
--  p = for i in 0..(rDim-1) list ((random QQ) + 13/29);
--  myJac = evalJac(p,phiCoords);
--  );
--  if (rank myJac == rDim) then print param;
--)

-- the output of the previous loop, saved in a list.
minimalProblems = (
(1, 0, 0, 1, 2, 0, 2, 2, 1),
(1, 0, 0, 2, 0, 0, 4, 0, 1),
(1, 0, 0, 2, 1, 0, 2, 1, 1),
(1, 0, 0, 3, 0, 0, 2, 0, 1),
(2, 0, 0, 0, 2, 0, 3, 1, 1),
(2, 0, 0, 0, 2, 0, 3, 2, 1),
(2, 0, 0, 1, 0, 1, 3, 1, 1),
(2, 1, 0, 0, 1, 0, 2, 1, 1),
(2, 1, 0, 1, 0, 0, 2, 0, 1),
(3, 0, 0, 0, 0, 2, 2, 2, 1),
(3, 0, 0, 0, 1, 0, 4, 1, 1),
(3, 0, 0, 0, 1, 1, 2, 1, 1),
(3, 0, 0, 0, 1, 1, 2, 2, 1),
(3, 0, 0, 0, 2, 0, 2, 1, 1),
(3, 0, 0, 0, 2, 0, 2, 2, 1),
(3, 0, 0, 1, 0, 0, 3, 0, 1),
(3, 0, 0, 1, 0, 1, 2, 1, 1),
(3, 0, 0, 1, 1, 0, 2, 1, 1),
(3, 0, 0, 2, 0, 0, 2, 0, 1),
(3, 1, 0, 0, 0, 0, 3, 0, 1),
(3, 2, 0, 0, 0, 0, 2, 0, 0),
(4, 0, 0, 0, 0, 0, 5, 0, 1),
(4, 1, 0, 0, 0, 0, 3, 0, 0),
(4, 1, 0, 0, 0, 0, 2, 0, 1),
(5, 0, 0, 0, 0, 0, 4, 0, 0),
(5, 0, 0, 0, 0, 1, 2, 1, 1),
(5, 0, 0, 0, 1, 0, 2, 1, 1),
(5, 0, 0, 1, 0, 0, 2, 0, 1),
(6, 1, 0, 0, 0, 0, 2, 0, 0),
(7, 0, 0, 0, 0, 0, 2, 0, 1),
(9, 0, 0, 0, 0, 0, 2, 0, 0)
)

--from now on, the code generates Julia scripts that run monodromy on 
--equations that have been computed in M2

--function that eats a minimal problem in form of a signature such
--as (1,0,0,3,0,0,2,0,1), and outputs a pair (Ringy, Leqs), where:
--Leqs is a list of fractions that are the coordinate functions of the
--rational joint camera map on an affine chart
--Ringy is the ring of the numerators/denoms of the fractions in Leqs 

setUpEquations = (a,b,c,d,e,f,m,incid,extra) -> (
  Leqs = mapCoords(a,b,c,d,e,f,m,incid,extra);
  Ringy = ring numerator Leqs_0;
  (Ringy, Leqs)
)

--helping function to remove the first and the last letter in a string
removeBrackets = str -> substring(1, (length str)-2, str)

--computes the numerator of a rational function, but also allows
--polynomials as input
num = F -> (
    if class class F === PolynomialRing then return F;
    numerator F
)

--computes the denominator of a rational function, but also allows
--polynomials as input
den = F -> (
    if class class F === PolynomialRing then return 1;
    denominator F
)

--input: a minimal problem, e.g. in form of a signature
--output: a string that is a Julia script to compute the degree of the
--joint camera map using monodromy, including interval certification.
setUpMonodromy = problem -> (
    (R,Leqs) = setUpEquations (problem);
    S = QQ[v_1..v_(#gens R),u_1..u_(#Leqs)];
    f = map(S,R,{v_1..v_(#gens R)});
    eqs = apply(#Leqs, i -> (f num Leqs#i) - u_(i+1)*(f den Leqs#i));
    
    str = "using HomotopyContinuation\n";
    variablesAndParams = apply(gens S, toString);
    str = str|"@var " | concatenate(mingle(variablesAndParams,
      #variablesAndParams:1)) | "\n";
    eqStr = toString eqs;
    variables = toString (v_1..v_(#gens R));
    params = toString (u_1..u_(#Leqs));
    str = str|"f = [" | removeBrackets(eqStr) | "]\n";
    str = str|"V = [" | removeBrackets(variables) | "];\n";
    str = str|"U = [" | removeBrackets(params) | "];\n";
    str = str|"F = System(f; variables = V, parameters = U)\n";
    
    str = str|"V0 = randn(ComplexF64, "|#Leqs|")\n";
    str = str|"F0 = System(f, variables = U, parameters = V)\n";
    str = str|"S0 = solve(F0, start_system = :total_degree,
      target_parameters = V0)\n";
    str = str|"U0 = solutions(S0)[1]\n";
    
    str = str|"R = monodromy_solve(F, V0, U0)\n";
    str = str|"cert = certify(F,R)\n";
    str
)

--helping function that determines the name of the file where the Julia
--script is stored.
signature = problem -> (
    str = "";
    scan(#problem, i -> str = str|(problem#i));
    str
)

--input: a minimal problem, e.g. in form of a signature
--this function outputs a Julia file with the string returned by
--setUpMonodromy of that problem.

-- comment out if directory already present
makeDirectory("./Minimal Problems/")

makeJuliaFile = problem -> (
    name = "./Minimal Problems/J"|(signature problem)|".jl";
    file = openOut name;
    file << setUpMonodromy(problem);
    file << endl;
    file << close;
)

-- Example:
-- makeJuliaFile (2,1,0,0,1,2,0,0,1)

--the resulting file can be run on in Julia like this:
--include("J100120201.jl")
--assuming the Julia interactive interface is run from the folder
--"Minimal Problems".

for problem in minimalProblems do (
  makeJuliaFile(problem);
)

