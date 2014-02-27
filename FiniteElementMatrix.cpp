#include "FiniteElementMatrix.h"
#include "bracket.h"

#include <cstdlib>
#include <ctime>


FiniteElementMatrix::FiniteElementMatrix(unsigned p,double simplex_peaks[4][3], double Eps[3][3], double Mu[3][3]) : m_P(p)
{
    for (unsigned i=0;i<4;i++)
    {
        for (unsigned j=0;j<3;j++)
        {
            m_Peaks[i][j]=simplex_peaks[i][j];
            if (i<=2)
            {
                m_MatrixMu[i][j]=Mu[i][j];
                m_MatrixEps[i][j]=Eps[i][j];
            }
        }
    }
}

FiniteElementMatrix::~FiniteElementMatrix()
{

}

Power_t FiniteElementMatrix::DefPowers(unsigned n,unsigned m)
{
    Power_t powers={0,0,0,0};

    if (n==1)
    {
        powers.p1=powers.p1+1;
    }
    if (n==2)
    {
        powers.p2=powers.p2+1;
    }
    if (n==3)
    {
        powers.p3=powers.p3+1;
    }
    if (n==4)
    {
        powers.p4=powers.p4+1;
    }

    if (m==1)
    {
        powers.p1=powers.p1+1;
    }
    if (m==2)
    {
        powers.p2=powers.p2+1;
    }
    if (m==3)
    {
        powers.p3=powers.p3+1;
    }
    if (m==4)
    {
        powers.p4=powers.p4+1;
    }
    return powers;
}

void FiniteElementMatrix::AddVTVProduct(std::vector<double> a, std::vector<double> b, std::vector<double> c,
                                        std::vector<double> d, double M[][3], unsigned n1,
                                        unsigned m1, unsigned n2, unsigned m2)
{
       double s_a, s_b;
       std::vector<double> M_a;
       std::vector<double> M_b;
       for (unsigned i=0;i<3;i++)
       {
           s_a=0; s_b=0;
           for (unsigned j=0;j<3;j++)
           {
               s_a=s_a+a[j]*M[j][i];
               s_b=s_b+b[j]*M[j][i];
           }

           M_a.push_back(s_a);
           M_b.push_back(s_b);
       }

       //������� ��������� ������ �������������

       double s;
       std::vector<double> gains;

       s=0;
       for (unsigned i=0;i<3;i++)
       {
           s=s+c[i]*M_a[i];
       }
       gains.push_back(s);

       s=0;
       for (unsigned i=0;i<3;i++)
       {
           s=s-d[i]*M_a[i];
       }
       gains.push_back(s);

       s=0;
       for (unsigned i=0;i<3;i++)
       {
           s=s-c[i]*M_b[i];
       }
       gains.push_back(s);

       s=0;
       for (unsigned i=0;i<3;i++)
       {
           s=s+d[i]*M_b[i];
       }
       gains.push_back(s);

       //��������� ������ ��������

       Power_t local_powers={0,0,0,0};
       std::vector<Power_t> powers;

       local_powers=DefPowers(n1,n2);

       powers.push_back(local_powers);

       //���������� ��� ����� n1,m2; m1,n2; m2,m2

       local_powers=DefPowers(n1,m2);

       powers.push_back(local_powers);

       local_powers=DefPowers(n2,m1);

       powers.push_back(local_powers);

       local_powers=DefPowers(m1,m2);

       powers.push_back(local_powers);

       //���������� ������ ��������� � ������
       AddToVectBracket(gains, powers);
}


void FiniteElementMatrix::MatrixInit(double eps[3][3], double mu[3][3])
{
    //Bracket b;
    unsigned ilow, ihigh,
        jlow, jhigh,
        klow, khigh,
        llow, lhigh;

    unsigned in_ilow, in_ihigh,
             in_jlow, in_jhigh,
             in_klow, in_khigh,
             in_llow, in_lhigh;

    for (unsigned gamma=1;gamma<=4;gamma++)
    {
        for (unsigned beta = gamma;beta<=4;beta++)
        {
            if ((gamma==1)||(beta==1))
            {
                ilow=0;
                ihigh=m_P;
            }
            else
            {
                ilow=1;
                ihigh=m_P+1;
            }
            if ((gamma==2)||(beta==2))
            {
                jlow=0;
                jhigh=m_P;
            }
            else
            {
                jlow=1;
                jhigh=m_P+1;
            }
            if ((gamma==3)||(beta==3))
            {
                klow=0;
                khigh=m_P;
            }
            else
            {
                klow=1;
                khigh=m_P+1;
            }
            if ((gamma==4)||(beta==4))
            {
                llow=0;
                lhigh=m_P;
            }
            else
            {
                llow=1;
                lhigh=m_P+1;
            }

            for (unsigned i=ilow;i<=ihigh;i++)
            {
                for (unsigned j=jlow;j<=jhigh;j++)
                {
                    for (unsigned k=klow;k<=khigh;k++)
                    {
                        for (unsigned l=llow;l<lhigh;l++)
                        {
                            if (i+j+k+l==m_P+2)
                            {
                                //����� ����� ������������� ���������� ����������;
                                //�������, ��� ���� ����� ������� ��� ��� ��������� �����

                                //� ������ ����� ����������� 4 ���������� ����������

                                //������ �������� ��������� �����

                                for (unsigned in_gamma=1;in_gamma<=4;in_gamma++)
                                {
                                    for (unsigned in_beta = in_gamma;in_beta<=4;in_beta++)
                                    {
                                        if ((in_gamma==1)||(in_beta==1))
                                        {
                                            in_ilow=0;
                                            in_ihigh=m_P;
                                        }
                                        else
                                        {
                                            in_ilow=1;
                                            in_ihigh=m_P+1;
                                        }
                                        if ((in_gamma==2)||(in_beta==2))
                                        {
                                            in_jlow=0;
                                            in_jhigh=m_P;
                                        }
                                        else
                                        {
                                            in_jlow=1;
                                            in_jhigh=m_P+1;
                                        }
                                        if ((in_gamma==3)||(in_beta==3))
                                        {
                                            in_klow=0;
                                            in_khigh=m_P;
                                        }
                                        else
                                        {
                                            in_klow=1;
                                            in_khigh=m_P+1;
                                        }
                                        if ((in_gamma==4)||(in_beta==4))
                                        {
                                            in_llow=0;
                                            in_lhigh=m_P;
                                        }
                                        else
                                        {
                                            in_llow=1;
                                            in_lhigh=m_P+1;
                                        }

                                        for (unsigned in_i=in_ilow;in_i<=in_ihigh;in_i++)
                                        {
                                            for (unsigned in_j=in_jlow;in_j<=in_jhigh;in_j++)
                                            {
                                                for (unsigned in_k=in_klow;in_k<=in_khigh;in_k++)
                                                {
                                                    for (unsigned in_l=in_llow;in_l<in_lhigh;in_l++)
                                                    {
                                                        if (in_i+in_j+in_k+in_l==m_P+2)
                                                        {
                                                            //����� ����� ������������� ���������� ����������;
                                                            //�������, ��� ���� ����� ������� ��� ��� ��������� �����

                                                            //� ������ ����� ����������� 4 ���������� ����������
                                                            //--------------------------------------------------------
                                                            //b.SetGains(gains);
                                                            //b.SetPowers(powers);
                                                            //--------------------------------------------------------

                                                            //������ ����� �������� � ������ ������� �������
                                                            unsigned n1, m1, n2, m2;
                                                            std::vector<double> a, b, c, d;
                                                            std::vector<unsigned> nm1, nm2;
                                                            nm1=Def_nm(gamma,beta);
                                                            nm2=Def_nm(in_gamma,in_beta);
                                                            n1=nm1[0];
                                                            m1=nm1[1];
                                                            n2=nm2[0];
                                                            m2=nm2[1];
                                                            a=DefVector(n1);
                                                            b=DefVector(m1);
                                                            c=DefVector(n2);
                                                            d=DefVector(m2);
                                                            AddVTVProduct(a, b, c, d, m_MatrixEps, n1, m1, n2, m2);
                                                        }
                                                    }//l ����������
                                                }//k ����������
                                            }//j ����������
                                        }//i ����������
                                    }//beta ����������
                                }//gamma ����������
                            }//�������� ������� (i+j+k+l==p+2)

                        }//l
                    }//k
                }//j
            }//i
        }//beta
    }//gamma
}//MatrixInit

void FiniteElementMatrix::AddToVectBracket(std::vector<double>& gains, std::vector<Power_t>& powers)
{
    m_VectBracket.push_back(Bracket(gains,powers));
}

void FiniteElementMatrix::ShowVectBracket()
{
    for (unsigned i=0;i<m_VectBracket.size();i++)
    {
        ((Bracket)(m_VectBracket.at(i))).ShowElements();
    }
}
//-----------------------����������� �������----------------------------------------------------------
std::vector<double> FiniteElementMatrix::DefVector(unsigned ind)
{
    std::vector<double> result;
    for (unsigned i=0;i<3;i++)
    {
        result.push_back(0);
    }

    std::vector<double> r1, r2, r3, r4;

    for (unsigned i=0;i<3;i++)
    {
        r1.push_back(m_Peaks[0][i]);
        r2.push_back(m_Peaks[1][i]);
        r3.push_back(m_Peaks[2][i]);
        r4.push_back(m_Peaks[3][i]);
    }
    std::vector<double> l1, l2, l3, l4;
    for (unsigned i=0;i<3;i++)
    {
        l1.push_back(r1[i]-r4[i]);
        l2.push_back(r2[i]-r4[i]);
        l3.push_back(r3[i]-r4[i]);
        l4.push_back(0);                   //????
    }
    double Icob;
    Icob=ScalarProduct(VectProduct(l1,l2),l3);


    if (ind==1)
    {
        result=VectProduct(l2,l3);
    }
    if (ind==2)
    {
        return VectProduct(l3,l4);
    }
    if (ind==3)
    {
        return VectProduct(l1,l2);
    }
    if (ind==4)
    {
        std::vector<double> a,b,c;
        a=VectProduct(l2,l3);
        b=VectProduct(l3,l4);
        c=VectProduct(l1,l2);
        for (unsigned i=0;i<3;i++)
        {
            result[i]=-a[i]-b[i]-c[i];
        }
    }
    for (unsigned i=0;i<3;i++)
    {
        result[i]=result[i]/Icob;
    }
    return result;
}
//----------------------------------------------------------------------------------------------------

//-----------------����������� �������� n, m ���� ������ gamma � beta---------------------------------
std::vector<unsigned> FiniteElementMatrix::Def_nm(unsigned gamma, unsigned beta)
{
    unsigned n, m;
    std::vector<unsigned> vect;
    if ((gamma==1)&&(beta==2))
    {
        n=3;
        m=4;
    }
    if ((gamma==1)&&(beta==3))
    {
        n=2;
        m=4;
    }
    if ((gamma==1)&&(beta==4))
    {
        n=2;
        m=3;
    }
    if ((gamma==2)&&(beta==3))
    {
        n=1;
        m=4;
    }
    if ((gamma==3)&&(beta==4))
    {
        n=1;
        m=2;
    }
    vect.push_back(n);
    vect.push_back(m);
    return vect;
}
//���������� ������, vect[0]=n, vect[1]=m; n<m;
//----------------------------------------------------------------------------------------------------
//----------------��������� ������������ ���� ��������------------------------------------------------
std::vector<double> FiniteElementMatrix::VectProduct(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> result;
    result.push_back(a[1]*b[2]-a[2]*b[1]);
    result.push_back(a[2]*b[0]-a[0]*b[2]);
    result.push_back(a[0]*b[1]-a[1]*b[0]);
    return result;
}
//----------------------------------------------------------------------------------------------------
//-------------��������� ������������ ���� ��������---------------------------------------------------
double FiniteElementMatrix::ScalarProduct(std::vector<double> a, std::vector<double> b)
{
    double result;
    result=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    return result;
}
//----------------------------------------------------------------------------------------------------
