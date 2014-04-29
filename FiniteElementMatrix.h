#ifndef FINITEELEMENTMATRIX_H
#define FINITEELEMENTMATRIX_H

#include "bracket.h"

#include <vector>
#include <list>
#include <stdlib.h>

class FiniteElementMatrix
{
public:
    static const unsigned m_CountPeaks=4;
    static const unsigned m_Dim=3;
    static const unsigned m_MaxOrderFact=25;
	static const unsigned m_QuadOrder=15;            //Порядок квадратурных формул
	//static const unsigned m_QuadOrder=2;

    FiniteElementMatrix(unsigned p, double simplex_peaks[m_CountPeaks][m_Dim], double Eps[m_Dim][m_Dim],
                        double Mu[m_Dim][m_Dim]);
    virtual ~FiniteElementMatrix();
    void MatrixInit();

	void NumMatrixInit();    //Формирование метрической матрицы и матрицы Эйлера с использованием численного интегрирования

    void AddToVectBracket(std::vector<GainPower_t>& terms, std::vector<Bracket>& vect_bracket);

    void ShowVectBracket(std::vector<Bracket>& vect_bracket);

    std::vector<Bracket> RotorCalc(Bracket& br, unsigned n, unsigned m);

    void ShowMatrixs();

    std::vector<Bracket> FormVectEigFunc(Bracket& br, unsigned n, unsigned m);

	void NumFormVectEigFunc(std::vector<Bracket>& vect_br, unsigned n, unsigned m, 
							double ksi1, double ksi2, double ksi3, 
							std::vector<double>& result);

	//void NumIntegration(std::vector<Bracket>& vect_br, std::vector<double>& result);

	double Integrate(Bracket& br);

private:

        unsigned m_P;

        std::vector<double> DefVector(unsigned ind);
        void Def_nm(unsigned gamma, unsigned beta, std::vector<unsigned>& vect);
        double m_Peaks[m_CountPeaks][m_Dim];       //массив координат (4) вершин тетраэдра

        double m_MatrixMu[m_Dim][m_Dim];    //Матрица тензора магнитной проницаемости
        double m_MatrixEps[m_Dim][m_Dim];    //Матрица тензора диэлектрической проницаемости

        void AddSilvester(unsigned gamma, unsigned beta, unsigned numb, unsigned ind, std::vector<Bracket>& vect_bracket);


        void findIndex(unsigned gamma, unsigned beta, unsigned order, unsigned& idxLow, unsigned& idxHigh);

        double *m_MetrMatrix;                    //Метрическая матрица
        double *m_EulerMatrix;                   //Матрица Эйлера
		double *m_NumEulerMatrix;
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

		//-----------------------------------------------------------------------------------------------------------------
		//Массивы корней полиномов Лежандра и весов для численного интегрирования
		std::vector<double> m_Weights;
		std::vector<double> m_Roots;

		void CompareMatrixs();
};

#endif // FINITEELEMENTMATRIX_H
