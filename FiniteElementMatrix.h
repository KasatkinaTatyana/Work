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
        FiniteElementMatrix(unsigned p, double simplex_peaks[4][3], double Eps[3][3], double Mu[3][3]);
        virtual ~FiniteElementMatrix();
        void MatrixInit(double eps[3][3], double mu[3][3]);

        std::vector<Bracket>& GetVectBracket() {return m_VectBracket;}

        void AddToVectBracket(std::vector<double>& gains, std::vector<Power_t>& powers);

        void ShowVectBracket();

        void AddVTVProduct(std::vector<double> a, std::vector<double> b, std::vector<double> c,
                           std::vector<double> d, double M[][3], unsigned n1, unsigned m1, unsigned n2, unsigned m2);
private:
        unsigned m_P;
        Power_t DefPowers(unsigned n,unsigned m);
        std::vector<Bracket> m_VectBracket;
        //void AddVTVProduct(std::vector<double> a, std::vector<double> b, std::vector<double> c,
        //                   std::vector<double> d, double M[3][3], unsigned n1, unsigned m1, unsigned n2, unsigned m2);
        std::vector<double> DefVector(unsigned ind);
        std::vector<unsigned> Def_nm(unsigned gamma, unsigned beta);
        double m_Peaks[4][3];       //массив координат (4) вершин тетраэдра

        std::vector<double> VectProduct(std::vector<double> a, std::vector<double> b);
        double ScalarProduct(std::vector<double> a, std::vector<double> b);

        double m_MatrixMu[3][3];    //Матрица тензора магнитной проницаемости
        double m_MatrixEps[3][3];    //Матрица тензора диэлектрической проницаемости

};
