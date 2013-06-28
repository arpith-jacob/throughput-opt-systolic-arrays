#ifndef _Lattice_h_
#define _Lattice_h_

#if defined(__cplusplus)
extern "C" {
#endif

extern void AffineHermite ( Lattice *A, Lattice **H, Matrix **U );
extern void AffineSmith ( Lattice *A, Lattice **U, Lattice **V, Lattice
                          **Diag );
extern Lattice *ChangeLatticeDimension ( Lattice *A, int dimension );
extern Lattice *EmptyLattice ( int dimension );
extern Lattice *ExtractLinearPart ( Lattice *A );
extern int FindHermiteBasisofDomain ( Polyhedron *A, Matrix **B );
extern Lattice *Homogenise ( Lattice *A, Bool Forward );
extern int intcompare ( const void *a, const void *b );
extern Bool isEmptyLattice ( Lattice *A );
extern Bool isfulldim ( Matrix *m );
extern Bool IsLattice ( Matrix *m );
extern Bool isLinear ( Lattice *A );
extern LatticeUnion *LatticeDifference ( Lattice *A, Lattice *B );
extern Lattice *LatticeImage ( Lattice *A, Matrix *M );
extern Bool LatticeIncludes ( Lattice *A, Lattice *B );
extern Lattice *LatticeIntersection ( Lattice *X, Lattice *Y );
extern Lattice *LatticePreimage ( Lattice *L, Matrix *G );
extern LatticeUnion *LatticeSimplify ( LatticeUnion *latlist );
extern LatticeUnion *LatticeUnion_Alloc ( void );
extern void LatticeUnion_Free ( LatticeUnion *Head );
extern void PrintLatticeUnion ( FILE *fp, char *format, LatticeUnion
                                *Head );
extern Bool sameAffinepart ( Lattice *A, Lattice *B );
extern Bool sameLattice ( Lattice *A, Lattice *B );
extern LatticeUnion *Lattice2LatticeUnion(Lattice *X,Lattice *Y);

#if defined(__cplusplus)
}
#endif

#endif /* _Lattice_h_ */
