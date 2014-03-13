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
    static const unsigned m_CountPeaks=4;
    static const unsigned m_Dim=3;
    static const unsigned m_MaxOrderFact=25;

    FiniteElementMatrix(unsigned p, double simplex_peaks[m_CountPeaks][m_Dim], double Eps[m_Dim][m_Dim],
                        double Mu[m_Dim][m_Dim]);
    virtual ~FiniteElementMatrix();
    void MatrixInit();

    std::vector<Bracket>& GetVectBracket() {return m_VectBracket;}

    void AddToVectBracket(std::vector<double>& gains, std::vector<Power_t>& powers, std::vector<Bracket>& vect_bracket);

    void ShowVectBracket(std::vector<Bracket>& vect_bracket);

    void AddRotMatrixProduct(Bracket vect1[m_Dim], Bracket vect2[m_Dim], double M[][m_Dim]);

private:

        unsigned m_P;

        Power_t DefPowers(unsigned n,unsigned m);
        std::vector<Bracket> m_VectBracket;       //эта скобка соответствует паре двух элементов
        std::vector<Bracket> m_CurVectBracket;    //эта скобка соответсвует текущему элементу и будет храниться
                                                  //при проходе по всем остальным элементам

        void AddVTVProduct(std::vector<double> a, std::vector<double> b, std::vector<double> c,
                           std::vector<double> d, double M[][m_Dim], unsigned n1, unsigned m1, unsigned n2, unsigned m2,
                           std::vector<Bracket>& vect_bracket);
        std::vector<double> DefVector(unsigned ind);
        std::vector<unsigned> Def_nm(unsigned gamma, unsigned beta);
        double m_Peaks[m_CountPeaks][m_Dim];       //массив координат (4) вершин тетраэдра

        std::vector<double> VectProduct(std::vector<double> a, std::vector<double> b);
        double ScalarProduct(std::vector<double> a, std::vector<double> b);

        double m_MatrixMu[m_Dim][m_Dim];    //Матрица тензора магнитной проницаемости
        double m_MatrixEps[m_Dim][m_Dim];    //Матрица тензора диэлектрической проницаемости

        void AddSilvester(unsigned gamma, unsigned beta, unsigned numb, unsigned ind, std::vector<Bracket>& vect_bracket);

        void LocalPowersChange(Power_t& local_powers, unsigned ind,unsigned value);

        double Integrate (Bracket_t& br);

        double** m_MetrMatrix;                    //Метрическая матрица
        double** m_EulerMatrix;                   //Матрица Эйлера
        unsigned m_MatrixSize;    //размерность метрической матрицы и матрицы Эйлера

        //Вычисление факториалов--------------------------------------------------------------------------------------------
        double Fact(unsigned N);

        double m_ArrayFact[m_MaxOrderFact+1];

        inline double CalcFact(unsigned N);
        //------------------------------------------------------------------------------------------------------------------
};
