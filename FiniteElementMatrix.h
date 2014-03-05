#ifndef FINITEELEMENTMATRIX_H
#define FINITEELEMENTMATRIX_H

#endif // FINITEELEMENTMATRIX_H

#include <vector>
#include <list>
#include "bracket.h"
#include <stdlib.h>
using namespace std;

typedef class Bracket Bracket_t;

class FiniteElementMatrix
{
public:
        FiniteElementMatrix(unsigned p, double simplex_peaks[4][3], double Eps[3][3], double Mu[3][3]);
        virtual ~FiniteElementMatrix();
        void MatrixInit(double eps[3][3], double mu[3][3]);

        std::vector<Bracket>& GetVectBracket() {return m_VectBracket;}

        void AddToVectBracket(std::vector<double>& gains, std::vector<Power_t>& powers, std::vector<Bracket>& vect_bracket);

        void ShowVectBracket(std::vector<Bracket>& vect_bracket);

private:
        unsigned m_P;

        Power_t DefPowers(unsigned n,unsigned m);
        std::vector<Bracket> m_VectBracket;       //эта скобка соответствует паре двух элементов
        std::vector<Bracket> m_CurVectBracket;    //эта скобка соответсвует текущему элементу и будет храниться
                                                  //при проходе по всем остальным элементам

        void AddVTVProduct(std::vector<double> a, std::vector<double> b, std::vector<double> c,
                           std::vector<double> d, double M[3][3], unsigned n1, unsigned m1, unsigned n2, unsigned m2,
                           std::vector<Bracket>& vect_bracket);
        std::vector<double> DefVector(unsigned ind);
        std::vector<unsigned> Def_nm(unsigned gamma, unsigned beta);
        double m_Peaks[4][3];       //массив координат (4) вершин тетраэдра

        std::vector<double> VectProduct(std::vector<double> a, std::vector<double> b);
        double ScalarProduct(std::vector<double> a, std::vector<double> b);

        double m_MatrixMu[3][3];    //Матрица тензора магнитной проницаемости
        double m_MatrixEps[3][3];    //Матрица тензора диэлектрической проницаемости

        void AddSilvester(unsigned gamma, unsigned beta, unsigned numb, unsigned ind, std::vector<Bracket>& vect_bracket);

        void LocalPowersChange(Power_t& local_powers, unsigned ind,unsigned value);

        double Integrate (Bracket_t& br);

        double* m_MetrMatrix;                    //Метрическая матрица
        double* m_EulerMatrix;                   //Матрица Эйлера


        //Вычисление факториалов--------------------------------------------------------------------------------------------
        double Fact(unsigned N);

        double m_ArrayFact[26];

        inline double CalcFact(unsigned N);
        //------------------------------------------------------------------------------------------------------------------
};
