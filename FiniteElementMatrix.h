#ifndef FINITEELEMENTMATRIX_H
#define FINITEELEMENTMATRIX_H

#endif // FINITEELEMENTMATRIX_H

#include <vector>
#include "bracket.h"

class FiniteElementMatrix
{
public:
        FiniteElementMatrix(unsigned p);
        virtual ~FiniteElementMatrix();
        void MatrixInit(double eps[3][3], double mu[3][3]);
        void AddVTVProduct(std::vector<float> a, std::vector<float> b, std::vector<float> c,
                           std::vector<float> d, double M[3][3], unsigned n1, unsigned m1, unsigned n2, unsigned m2);;

private:
        unsigned m_P;
        Power_t FormPowers(unsigned n,unsigned m);

};
