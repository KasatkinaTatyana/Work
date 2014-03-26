#ifndef FINITEELEMENTMATRIX_H
#define FINITEELEMENTMATRIX_H

#endif // FINITEELEMENTMATRIX_H

#include <vector>
#include <list>
#include "bracket.h"
#include <stdlib.h>

/*struct Row_t
{
    std::vector<double> rows;
};*/


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

    void AddToVectBracket(std::vector<double>& gains, std::vector<Power_t>& powers, std::vector<Bracket>& vect_bracket);

    void ShowVectBracket(std::vector<Bracket>& vect_bracket);

    Bracket GeneralVectorTensorVectorProduct(std::vector<Bracket>& vect1, std::vector<Bracket>& vect2, double M[][m_Dim]);

    std::vector<Bracket> RotorCalc(Bracket& br, unsigned n, unsigned m);
    std::vector<double> ProductGradKsi(unsigned ind1, unsigned ind2);

    std::vector<Bracket> VectBracketProduct(std::vector<Bracket>& a, std::vector<Bracket>& b);

    void ShowMatrixs();

    std::vector<Bracket> FormVectEigFunc(Bracket& br, unsigned n, unsigned m);

private:

        unsigned m_P;

        Power_t DefPowers(unsigned n,unsigned m);
        std::vector<Bracket> m_InnerVectBracket;  //эта скобка соответствует второму элементу (внутренний цикл)
        std::vector<Bracket> m_CurVectBracket;    //эта скобка соответсвует текущему элементу и будет храниться
                                                  //при проходе по всем остальным элементам

        std::vector<double> DefVector(unsigned ind);
        std::vector<unsigned> Def_nm(unsigned gamma, unsigned beta);
        double m_Peaks[m_CountPeaks][m_Dim];       //массив координат (4) вершин тетраэдра

        std::vector<double> VectProduct(std::vector<double> a, std::vector<double> b);
        double ScalarProduct(std::vector<double> a, std::vector<double> b);

        double m_MatrixMu[m_Dim][m_Dim];    //Матрица тензора магнитной проницаемости
        double m_MatrixEps[m_Dim][m_Dim];    //Матрица тензора диэлектрической проницаемости

        void AddSilvester(unsigned gamma, unsigned beta, unsigned numb, unsigned ind, std::vector<Bracket>& vect_bracket);

        void LocalPowersChange(Power_t& local_powers, unsigned ind,unsigned value);

        double Integrate (Bracket& br);

        void findIndex(unsigned gamma, unsigned beta, unsigned order, unsigned& idxLow, unsigned& idxHigh);

        double *m_MetrMatrix;                    //Метрическая матрица
        double *m_EulerMatrix;                   //Матрица Эйлера
        unsigned m_MatrixSize;    //размерность метрической матрицы и матрицы Эйлера

        //Вычисление факториалов--------------------------------------------------------------------------------------------
        double Fact(unsigned N);

        double m_ArrayFact[m_MaxOrderFact+1];

        inline double CalcFact(unsigned N);
        //------------------------------------------------------------------------------------------------------------------
        //градиенты барицентрических координат
        std::vector<double> m_GradKsi_1;
        std::vector<double> m_GradKsi_2;
        std::vector<double> m_GradKsi_3;
        std::vector<double> m_GradKsi_4;

        //Значения унитарных векторов, направленных вдоль ребер тетраэдра
        std::vector<double> m_UnVect_1;
        std::vector<double> m_UnVect_2;
        std::vector<double> m_UnVect_3;

        //Якобиан замены координат
        double m_Icob;
};
