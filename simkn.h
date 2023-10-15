/* Copyright (C) 2023 blackedout01

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 
*/

#ifndef SIMKN_NO_TYPE_DEF
#include <stdint.h>
typedef double f64;
typedef uint64_t u64;
typedef uint32_t u32;
#endif

static void SkAssign3(f64 *LHS, f64 *RHS);

static f64 SkDot3(f64 *A, f64 *B);

static void SkAssignCross(f64 *LHS, f64 *A, f64 *B);

static f64 SkScalarTripleProduct(f64 *A, f64 *B, f64 *C);

static void SkDivide3(f64 *LHS, f64 Denom);

static void SkScale3(f64 *LHS, f64 S);

static void SkAdd3(f64 *LHS, f64 *Addend);

static f64 SkDeterminant3x3(f64 *E);

static void SkInverse3x3(f64 *E, f64 *O);

static f64 SkTet3InertiaMoment(f64 *P, u64 I);

static f64 SkTet3IntertiaProduct(f64 *P, u64 I, u64 J);

static void SkComputeInertia3x3(void *User, u64 TriangleCount, f64 Density, f64 *I0, f64 *CoM, f64 *M,
                                void (*GetTrianglePositions)(void *, u64, f64 *));

#ifdef SIMKN_IMPLEMENTATION

// NOTE(blackedout):
// Matrix elements are assumed to be stored in row major format by default.
#ifdef SIMKN_COL_MAJOR
#define X3(Row, Col) (3*(Col) + (Row))
#else
#define X3(Row, Col) (3*(Row) + (Col))
#endif

#define SkSq(X) ((X)*(X))

static void SkAssign3(f64 *LHS, f64 *RHS) {
    LHS[0] = RHS[0];
    LHS[1] = RHS[1];
    LHS[2] = RHS[2];
}

static f64 SkDot3(f64 *A, f64 *B) {
    f64 Result = A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    return Result;
}

static void SkAssignCross(f64 *LHS, f64 *A, f64 *B) {
    LHS[0] = A[1]*B[2] - B[1]*A[2];
    LHS[1] = A[2]*B[0] - B[2]*A[0];
    LHS[2] = A[0]*B[1] - B[0]*A[1];
}

static f64 SkScalarTripleProduct(f64 *A, f64 *B, f64 *C) {
    f64 Tmp[3];
    SkAssignCross(Tmp, B, C);
    f64 Result = SkDot3(A, Tmp);
    return Result;
}

static void SkScale3(f64 *LHS, f64 S) {
    LHS[0] *= S;
    LHS[1] *= S;
    LHS[2] *= S;
}

static void SkDivide3(f64 *LHS, f64 Denom) {
    LHS[0] /= Denom;
    LHS[1] /= Denom;
    LHS[2] /= Denom;
}

static void SkAdd3(f64 *LHS, f64 *Addend) {
    LHS[0] += Addend[0];
    LHS[1] += Addend[1];
    LHS[2] += Addend[2];
}

static f64 SkDeterminant3x3(f64 *E) {
    // NOTE(blackedout): Compare to SkCross and
    // https://matrixcalc.org/en/#determinant({{x_0,a_0,b_0},{x_1,a_1,b_1},{x_2,a_2,b_2}})
    
    f64 X0 = E[X3(0, 0)];
    f64 X1 = E[X3(1, 0)];
    f64 X2 = E[X3(2, 0)];
    
    f64 A0 = E[X3(0, 1)];
    f64 A1 = E[X3(1, 1)];
    f64 A2 = E[X3(2, 1)];
    
    f64 B0 = E[X3(0, 2)];
    f64 B1 = E[X3(1, 2)];
    f64 B2 = E[X3(2, 2)];
    
    f64 Result =
        X0*(A1*B2 - B1*A2) +
        X1*(A2*B0 - B2*A0) +
        X2*(A0*B1 - B0*A1);
    return Result;
}

static void SkInverse3x3(f64 *E, f64 *O) {
    // NOTE(blackedout): Compare to
    // https://matrixcalc.org/en/#inverse({{x_0,a_0,b_0},{x_1,a_1,b_1},{x_2,a_2,b_2}})
    
    f64 X0 = E[X3(0, 0)];
    f64 X1 = E[X3(1, 0)];
    f64 X2 = E[X3(2, 0)];
    
    f64 A0 = E[X3(0, 1)];
    f64 A1 = E[X3(1, 1)];
    f64 A2 = E[X3(2, 1)];
    
    f64 B0 = E[X3(0, 2)];
    f64 B1 = E[X3(1, 2)];
    f64 B2 = E[X3(2, 2)];
    
    f64 Det = SkDeterminant3x3(E);
    
    O[X3(0, 0)] = (A1*B2 - B1*A2)/Det;
    O[X3(0, 1)] = (A2*B0 - B2*A0)/Det;
    O[X3(0, 2)] = (A0*B1 - B0*A1)/Det;
    
    O[X3(1, 0)] = (B1*X2 - X1*B2)/Det;
    O[X3(1, 1)] = (B2*X0 - X2*B0)/Det;
    O[X3(1, 2)] = (B0*X1 - X0*B1)/Det;
    
    O[X3(2, 0)] = (X1*A2 - A1*X2)/Det;
    O[X3(2, 1)] = (X2*A0 - A2*X0)/Det;
    O[X3(2, 2)] = (X0*A1 - A0*X1)/Det;
}

static f64 SkTet3InertiaMoment(f64 *P, u64 I) {
    f64 Result = SkSq(P[3*0 + I]) + P[3*1 + I]*P[3*2 + I]
               + SkSq(P[3*1 + I]) + P[3*0 + I]*P[3*2 + I]
               + SkSq(P[3*2 + I]) + P[3*0 + I]*P[3*1 + I];
    return Result;
}

static f64 SkTet3IntertiaProduct(f64 *P, u64 I, u64 J) {
    f64 Result = 2.0*P[3*0 + I]*P[3*0 + J] + P[3*1 + I]*P[3*2 + J] + P[3*2 + I]*P[3*1 + J]
               + 2.0*P[3*1 + I]*P[3*1 + J] + P[3*0 + I]*P[3*2 + J] + P[3*2 + I]*P[3*0 + J]
               + 2.0*P[3*2 + I]*P[3*2 + J] + P[3*0 + I]*P[3*1 + J] + P[3*1 + I]*P[3*0 + J];
    return Result;
}

static void SkCuboidInertia3(f64 Density, f64 *Size, f64 *I) {
    f64 XX = SkSq(Size[0]);
    f64 YY = SkSq(Size[1]);
    f64 ZZ = SkSq(Size[2]);
    f64 Mass = Density*Size[0]*Size[1]*Size[2];
    I[0] = Mass*(YY + ZZ)/12.0;
    I[1] = Mass*(XX + ZZ)/12.0;
    I[2] = Mass*(XX + YY)/12.0;
}

static void SkComputeInertia3x3(void *User, u64 TriangleCount, f64 Density, f64 *I0, f64 *CoM, f64 *M,
                                void (*GetTrianglePositions)(void *, u64, f64 *)) {
    f64 Mass = 0.0;
    f64 MassCenter[3] = {0};
    f64 Ia = 0.0, Ib = 0.0, Ic = 0.0, Iap = 0.0, Ibp = 0.0, Icp = 0.0;
    f64 P[9];
    for(u64 I = 0; I < TriangleCount; ++I) {
        GetTrianglePositions(User, I, P);
        
        // NOTE(blackedout): The following three properties are signed
        f64 DetJ = SkScalarTripleProduct(P + 0, P + 3, P + 6);
        f64 TetVolume = DetJ/6.0;
        f64 TetMass = Density*TetVolume;
        
        f64 TetMassCenter[3] = {0};
        SkAdd3(TetMassCenter, P + 0);
        SkAdd3(TetMassCenter, P + 3);
        SkAdd3(TetMassCenter, P + 6);
        SkDivide3(TetMassCenter, 4.0);
        
        f64 V100 = SkTet3InertiaMoment(P, 0);
        f64 V010 = SkTet3InertiaMoment(P, 1);
        f64 V001 = SkTet3InertiaMoment(P, 2);
        
        Ia += DetJ*(V010 + V001);
        Ib += DetJ*(V100 + V001);
        Ic += DetJ*(V100 + V010);
        Iap += DetJ*SkTet3IntertiaProduct(P, 1, 2);
        Ibp += DetJ*SkTet3IntertiaProduct(P, 0, 1);
        Icp += DetJ*SkTet3IntertiaProduct(P, 0, 2);
        
        SkScale3(TetMassCenter, TetMass);
        SkAdd3(MassCenter, TetMassCenter);
        Mass += TetMass;
    }
    
    SkDivide3(MassCenter, Mass);
    Ia = Density*Ia/60.0 - Mass*(SkSq(MassCenter[1]) + SkSq(MassCenter[2]));
    Ib = Density*Ib/60.0 - Mass*(SkSq(MassCenter[0]) + SkSq(MassCenter[2]));
    Ic = Density*Ic/60.0 - Mass*(SkSq(MassCenter[0]) + SkSq(MassCenter[1]));
    Iap = Density*Iap/120.0 - Mass*(MassCenter[1]*MassCenter[2]);
    Ibp = Density*Ibp/120.0 - Mass*(MassCenter[0]*MassCenter[1]);
    Icp = Density*Icp/120.0 - Mass*(MassCenter[0]*MassCenter[2]);
    
    I0[0] = Ia;
    I0[4] = Ib;
    I0[8] = Ic;
    I0[1] = I0[3] = -Ibp;
    I0[2] = I0[6] = -Icp;
    I0[5] = I0[7] = -Iap;
    
    SkAssign3(CoM, MassCenter);
    *M = Mass;
}

#endif
