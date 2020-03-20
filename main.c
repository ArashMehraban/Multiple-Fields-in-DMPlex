const char help[] = "Test multiple fields in PETSc\n";

#include<petscdmplex.h>
#include "petscfe.h"          
#include "petscdmlabel.h"     
#include "petscds.h" 

PetscErrorCode PetscFECreateByDegree(DM dm, PetscInt dim, PetscInt Nc, PetscBool isSimplex, const char prefix[], PetscInt order, PetscBool continuous, PetscFE *fem); 

int main(int argc, char **argv) {
  PetscInt     ierr;
  char         meshFile[PETSC_MAX_PATH_LEN];
  PetscInt     degree;       //command line polynomial degree for displacement 
  PetscInt     degreeP = 0;  // degreeP = 0 create P_(-1)
  PetscBool    continuousDisplacement = PETSC_TRUE; // Displacement is continuous
  PetscBool    continuousPressure = PETSC_FALSE;    // Pressure is discontinuous
  DM           dm;
  DM           distributedMesh = NULL;
  PetscFE      fe, feP;
  PetscInt     ncompu = 3;
  PetscInt     ncompp = 4;
  PetscInt     dim = 3; 
  Vec          vec;
  PetscInt     vec_sz;
  PetscBool    interpolate = PETSC_FALSE; 
  
  ierr = PetscInitialize(&argc, &argv, NULL, help);
  if (ierr)
    return ierr;

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "Displacement/Pressure fields", NULL); CHKERRQ(ierr);
    ierr = PetscOptionsString("-mesh", "Read mesh from file", NULL, meshFile, meshFile, sizeof(meshFile), NULL);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-degree", "Polynomial degree of tensor product basis",NULL, degree, &degree,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);
   
  if(degree >1)
     interpolate = PETSC_TRUE;
  ierr = DMPlexCreateFromFile(PETSC_COMM_WORLD, meshFile, interpolate, &dm); CHKERRQ(ierr);
  ierr = DMPlexDistribute(dm, 0, NULL, &distributedMesh);
  if (distributedMesh) {
     DMDestroy(&dm);
     dm  = distributedMesh;
  }

  ierr = PetscFECreateByDegree(dm, dim, ncompu, PETSC_FALSE, NULL, degree, continuousDisplacement, &fe);
  ierr=  DMAddField(dm, NULL, (PetscObject)fe);

  ierr = PetscFECreateByDegree(dm, dim, ncompp, PETSC_FALSE, NULL, degreeP, continuousPressure, &feP);  
  ierr = DMAddField(dm, NULL, (PetscObject)feP); 

  ierr = DMCreateDS(dm);
  ierr = DMPlexSetClosurePermutationTensor(dm, PETSC_DETERMINE, NULL);
  ierr = DMView(dm,PETSC_VIEWER_STDOUT_WORLD);
  ierr = DMGetLocalVector(dm, &vec);
  ierr = VecGetSize(vec, &vec_sz);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Vec size: %D\n", vec_sz);
  
  

  //cleanup
  DMDestroy(&dm);
  PetscFEDestroy(&fe);
  PetscFEDestroy(&feP);

  return PetscFinalize();
}


PetscErrorCode PetscFECreateByDegree(DM dm, PetscInt dim, PetscInt Nc,
                                      PetscBool isSimplex, const char
prefix[],
                                      PetscInt order, PetscBool continuous,
                                      PetscFE *fem) {
   PetscQuadrature q, fq;
   DM              K;
   PetscSpace      P;
   PetscDualSpace  Q;
   PetscInt        quadPointsPerEdge;
   PetscBool       tensor = isSimplex ? PETSC_FALSE : PETSC_TRUE;
   PetscErrorCode  ierr;

   PetscFunctionBeginUser;
   /* Create space */
   ierr = PetscSpaceCreate(PetscObjectComm((PetscObject) dm), &P);
CHKERRQ(ierr);
   ierr = PetscObjectSetOptionsPrefix((PetscObject) P, prefix);
CHKERRQ(ierr);
   ierr = PetscSpacePolynomialSetTensor(P, tensor); CHKERRQ(ierr);
   ierr = PetscSpaceSetFromOptions(P); CHKERRQ(ierr);
   ierr = PetscSpaceSetNumComponents(P, Nc); CHKERRQ(ierr);
   ierr = PetscSpaceSetNumVariables(P, dim); CHKERRQ(ierr);
   ierr = PetscSpaceSetDegree(P, order, order); CHKERRQ(ierr);
   ierr = PetscSpaceSetUp(P); CHKERRQ(ierr);
   ierr = PetscSpacePolynomialGetTensor(P, &tensor); CHKERRQ(ierr);
   /* Create dual space */
   ierr = PetscDualSpaceCreate(PetscObjectComm((PetscObject) dm), &Q);
   CHKERRQ(ierr);
   ierr = PetscDualSpaceSetType(Q,PETSCDUALSPACELAGRANGE); CHKERRQ(ierr);
   if (!continuous) {
     ierr = PetscDualSpaceLagrangeSetContinuity(Q, PETSC_FALSE);
CHKERRQ(ierr);
   }
   ierr = PetscObjectSetOptionsPrefix((PetscObject) Q, prefix);
CHKERRQ(ierr);
   ierr = PetscDualSpaceCreateReferenceCell(Q, dim, isSimplex, &K);
CHKERRQ(ierr);
   ierr = PetscDualSpaceSetDM(Q, K); CHKERRQ(ierr);
   ierr = DMDestroy(&K); CHKERRQ(ierr);
   ierr = PetscDualSpaceSetNumComponents(Q, Nc); CHKERRQ(ierr);
   ierr = PetscDualSpaceSetOrder(Q, order); CHKERRQ(ierr);
   ierr = PetscDualSpaceLagrangeSetTensor(Q, tensor); CHKERRQ(ierr);
   ierr = PetscDualSpaceSetFromOptions(Q); CHKERRQ(ierr);
   ierr = PetscDualSpaceSetUp(Q); CHKERRQ(ierr);
   /* Create element */
   ierr = PetscFECreate(PetscObjectComm((PetscObject) dm), fem);
CHKERRQ(ierr);
   ierr = PetscObjectSetOptionsPrefix((PetscObject) *fem, prefix);
CHKERRQ(ierr);
   ierr = PetscFESetFromOptions(*fem); CHKERRQ(ierr);
   ierr = PetscFESetBasisSpace(*fem, P); CHKERRQ(ierr);
   ierr = PetscFESetDualSpace(*fem, Q); CHKERRQ(ierr);
   ierr = PetscFESetNumComponents(*fem, Nc); CHKERRQ(ierr);
   ierr = PetscFESetUp(*fem); CHKERRQ(ierr);
   ierr = PetscSpaceDestroy(&P); CHKERRQ(ierr);
   ierr = PetscDualSpaceDestroy(&Q); CHKERRQ(ierr);
   /* Create quadrature */
   quadPointsPerEdge = PetscMax(order + 1,1);
   if (isSimplex) {
     ierr = PetscDTStroudConicalQuadrature(dim,   1, quadPointsPerEdge,
-1.0, 1.0,
                                           &q); CHKERRQ(ierr);
     ierr = PetscDTStroudConicalQuadrature(dim-1, 1, quadPointsPerEdge,
-1.0, 1.0,
                                           &fq); CHKERRQ(ierr);
   } else {
     ierr = PetscDTGaussTensorQuadrature(dim,   1, quadPointsPerEdge,
-1.0, 1.0,
                                         &q); CHKERRQ(ierr);
     ierr = PetscDTGaussTensorQuadrature(dim-1, 1, quadPointsPerEdge,
-1.0, 1.0,
                                         &fq); CHKERRQ(ierr);
   }
   ierr = PetscFESetQuadrature(*fem, q); CHKERRQ(ierr);
   ierr = PetscFESetFaceQuadrature(*fem, fq); CHKERRQ(ierr);
   ierr = PetscQuadratureDestroy(&q); CHKERRQ(ierr);
   ierr = PetscQuadratureDestroy(&fq); CHKERRQ(ierr);

   PetscFunctionReturn(0);
}



