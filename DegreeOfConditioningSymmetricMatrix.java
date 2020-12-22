package vsu.labs;

import org.apache.commons.math3.linear.*;

import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.math3.linear.RealMatrix;

import static org.apache.commons.math3.linear.MatrixUtils.createRealMatrix;
import static org.apache.commons.math3.linear.MatrixUtils.inverse;

public class DegreeOfConditioningSymmetricMatrix
{
    static int n = 10;

    public static RealMatrix initL(int eps)
    {
        double[][] L1 = new double[n][n];
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                if (i == j)
                {
                    L1[i][j] = Math.random() * 2 * Math.pow(10, eps) - Math.pow(10, eps);;
                }
                else
                {
                    L1[i][j] = 0;
                }
            }
        }
        RealMatrix L = createRealMatrix(L1);//диагональная матрица с собственными значениями

        return L;
    }

    public static RealMatrix initH()
    {
        int w = (int) (Math.random() * n); //число [0; n)
        double[][] W1 = new double[n][1];
        for(int i = 0; i < n; i++)
        {
            if (i == w)
            {
                W1[i][0] = 1;
            }
            else
            {
                W1[i][0] = 0;
            }
        }
        RealMatrix W = createRealMatrix(W1);//случайный единичный вектор столбец W

        double[][] E1 = new double[n][n];
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                if (i == j)
                {
                    E1[i][j] = 1;
                }
                else
                {
                    E1[i][j] = 0;
                }
            }
        }
        RealMatrix E = createRealMatrix(E1);//единичная матрица

        RealMatrix Wt = W.transpose();
        RealMatrix WmultWt = W.multiply(Wt);
        RealMatrix H = E.subtract(WmultWt.scalarMultiply(2));//Матрица Хаусхолдера: H = E - 2 * W * Wt
        //столбцы этой матрицы - соответствующие собственные векторы

        return H;
    }

    public static RealMatrix initA(RealMatrix L, RealMatrix H)
    {
        RealMatrix Ht = H.transpose();
        RealMatrix HmultL = H.multiply(L);
        RealMatrix A = HmultL.multiply(Ht);//Матрица A = H * L * Ht
        return A;
    }

    public static double directIterationMethod(RealMatrix A) //Метод прямых итераций
    {
        double eps = 0.001;

        double[][] xCurr1 = new double[n][1];

        for(int i = 0; i < n; i++)
        {
            xCurr1[i][0] = 1;
        }
        RealMatrix xCurr = createRealMatrix(xCurr1);//произвольный вектор (в этом случае состоит из единиц)

        int k = 0;

        RealMatrix v = xCurr.scalarMultiply(1/xCurr.getNorm());//v[k] = x[k]/||x[k]||
        RealMatrix xNext = A.multiply(v);//x[k+1] = A * v[k]
        RealMatrix vt = v.transpose();
        RealMatrix sCurr = vt.multiply(xNext); //s[k] = vt[k] * x[k+1]

        xCurr = xNext.copy();
        RealMatrix sPrev = sCurr.copy();

        boolean flag = true;
        double Ln = 0;
        RealMatrix xn = createRealMatrix(new double[n][1]);

        while (flag)
        {
            k++;
            v = xCurr.scalarMultiply(1/xCurr.getNorm());
            xNext = A.multiply(v);
            vt = v.transpose();
            sCurr = vt.multiply(xNext);

            if(Math.abs(sCurr.getEntry(0, 0) - sPrev.getEntry(0, 0)) <= eps)//|s[k] - s[k-1]| <= eps
            {
                Ln = sCurr.getNorm(); // Ln = s[k] - максимальное по модулю собственное значение
                xn = v.copy(); // xn = v[k] - соответствующий собственный вектор
                flag = false;
            }

            xCurr = xNext.copy();
            sPrev = sCurr.copy();
        }

        System.out.println(Ln);
        System.out.println(xn);
        System.out.println(k);

        return Ln;
    }

    public static double reverseIterationMethod(RealMatrix A) //Метод обратных итераций
    {
        double eps = 0.001;

        double[][] xCurr1 = new double[n][1];

        for(int i = 0; i < n; i++)
        {
            xCurr1[i][0] = 1;
        }
        RealMatrix xCurr = createRealMatrix(xCurr1);//произвольный вектор (в этом случае состоит из единиц)

        int k = 0;

        RealMatrix v = xCurr.scalarMultiply(1/xCurr.getNorm());//v[k] = x[k]/||x[k]||
        RealMatrix xNext = (MatrixUtils.inverse(A)).multiply(v);//x[k+1] = A[-1](обратная матрица) * v[k]
        RealMatrix vt = v.transpose();
        RealMatrix aCurr = vt.multiply(xNext); //a[k] = vt[k] * x[k+1]

        xCurr = xNext.copy();
        RealMatrix aPrev = aCurr.copy();

        boolean flag = true;
        double L1 = 0;
        RealMatrix x1 = createRealMatrix(new double[n][1]);

        while (flag)
        {
            k++;
            v = xCurr.scalarMultiply(1/xCurr.getNorm());
            xNext = (MatrixUtils.inverse(A)).multiply(v);
            vt = v.transpose();
            aCurr = vt.multiply(xNext);

            if(Math.abs(aCurr.getEntry(0, 0) - aPrev.getEntry(0, 0)) <= eps)//|a[k] - a[k-1]| <= eps
            {
                L1 = 1 / aCurr.getNorm(); // L1 = 1 / a[k] - максимальное по модулю собственное значение обратной матрицы
                x1 = v.copy(); // x1 = v[k] - соответствующий собственный вектор
                flag = false;
            }

            xCurr = xNext.copy();
            aPrev = aCurr.copy();
        }

        System.out.println(L1);
        System.out.println(x1);
        System.out.println(k);

        return 1 / L1;
    }

    public static double conditionNumberOfTheMatrix(double L1, double Ln)
    {
        return Math.abs(Ln / L1);
    }



}
