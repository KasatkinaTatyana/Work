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

    std::vector<Bracket> RotorCalc(Bracket& br, unsigned n, unsigned m);

    void ShowMatrixs();

    std::vector<Bracket> FormVectEigFunc(Bracket& br, unsigned n, unsigned m);

private:

        unsigned m_P;

        Power_t DefPowers(unsigned n,unsigned m);
        std::vector<Bracket> m_InnerVectBracket;  //��� ������ ������������� ������� �������� (���������� ����)
        std::vector<Bracket> m_CurVectBracket;    //��� ������ ������������ �������� �������� � ����� ���������
                                                  //��� ������� �� ���� ��������� ���������

        std::vector<double> DefVector(unsigned ind);
        std::vector<unsigned> Def_nm(unsigned gamma, unsigned beta);
        double m_Peaks[m_CountPeaks][m_Dim];       //������ ��������� (4) ������ ���������

        double m_MatrixMu[m_Dim][m_Dim];    //������� ������� ��������� �������������
        double m_MatrixEps[m_Dim][m_Dim];    //������� ������� ��������������� �������������

        void AddSilvester(unsigned gamma, unsigned beta, unsigned numb, unsigned ind, std::vector<Bracket>& vect_bracket);

        void LocalPowersChange(Power_t& local_powers, unsigned ind,unsigned value);

        double Integrate (Bracket& br);

        void findIndex(unsigned gamma, unsigned beta, unsigned order, unsigned& idxLow, unsigned& idxHigh);

        double *m_MetrMatrix;                    //����������� �������
        double *m_EulerMatrix;                   //������� ������
        unsigned m_MatrixSize;    //����������� ����������� ������� � ������� ������

        //���������� �����������--------------------------------------------------------------------------------------------
        double Fact(unsigned N);

        double m_ArrayFact[m_MaxOrderFact+1];

        inline double CalcFact(unsigned N);
        //------------------------------------------------------------------------------------------------------------------
        //��������� ���������������� ���������
        std::vector<double> m_GradKsi_1;
        std::vector<double> m_GradKsi_2;
        std::vector<double> m_GradKsi_3;
        std::vector<double> m_GradKsi_4;

        //�������� ��������� ��������, ������������ ����� ����� ���������
        std::vector<double> m_UnVect_1;
        std::vector<double> m_UnVect_2;
        std::vector<double> m_UnVect_3;

        //������� ������ ���������
        double m_Icob;
};
