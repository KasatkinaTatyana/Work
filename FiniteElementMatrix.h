#ifndef FINITEELEMENTMATRIX_H
#define FINITEELEMENTMATRIX_H

#endif // FINITEELEMENTMATRIX_H

#include <vector>
#include <list>
#include "bracket.h"
#include <stdlib.h>
using namespace std;

class FiniteElementMatrix
{
public:
        FiniteElementMatrix(unsigned p);
        virtual ~FiniteElementMatrix();
        void MatrixInit(double eps[3][3], double mu[3][3]);
        void AddVTVProduct(std::vector<double> a, std::vector<double> b, std::vector<double> c,
                           std::vector<double> d, double M[3][3], unsigned n1, unsigned m1, unsigned n2, unsigned m2);

        std::vector<Bracket>& GetVectBracket() {return m_VectBracket;}

        void AddToVectBracket(std::vector<double>& gains, std::vector<Power_t>& powers);

        void ShowVectBracket();

private:
        unsigned m_P;
        Power_t FormPowers(unsigned n,unsigned m);
        std::vector<Bracket> m_VectBracket;

};
