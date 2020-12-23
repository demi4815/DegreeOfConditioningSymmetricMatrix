package vsu.labs;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.Precision;

public class Tests
{
    public static void test1()
    {
        int eps  = 2;
        RealMatrix L = DegreeOfConditioningSymmetricMatrix.initL(eps);
        RealMatrix H = DegreeOfConditioningSymmetricMatrix.initH();
        RealMatrix A = DegreeOfConditioningSymmetricMatrix.initA(L, H);

        System.out.println(L);
        System.out.println(H);
        System.out.println(A);

        System.out.println();

        double Lmax = DegreeOfConditioningSymmetricMatrix.getMaxL(L);
        double Lmin = DegreeOfConditioningSymmetricMatrix.getMinL(L, eps);
        double m1 = DegreeOfConditioningSymmetricMatrix.conditionNumberOfTheMatrix(Lmin, Lmax);

        System.out.println(Lmax);
        System.out.println(Lmin);
        System.out.println(m1);
        System.out.println();

        double Ln = DegreeOfConditioningSymmetricMatrix.directIterationMethod(A);
        System.out.println();
        double L1 = DegreeOfConditioningSymmetricMatrix.reverseIterationMethod(A);

        System.out.println();

        System.out.println(Ln);
        System.out.println(L1);
        double m2 = DegreeOfConditioningSymmetricMatrix.conditionNumberOfTheMatrix(L1, Ln);
        System.out.println(m2);
    }

    public static void test2()
    {
        int IER = 1;
        double maxN = 30;
        for (int n = 10; n <= maxN; n += 10)
        {
            for (int eps = 1; eps <= 2; eps += 1)
            {
                double absL1 = 0;
                double absLn = 0;
                double absM = 0;

                double otnosL1 = 0;
                double otnosLn = 0;
                double otnosM = 0;

                for (int cnt = 0; cnt < 10; cnt++)
                {
                    DegreeOfConditioningSymmetricMatrix.n = n;

                    RealMatrix L = DegreeOfConditioningSymmetricMatrix.initL(eps);
                    RealMatrix H = DegreeOfConditioningSymmetricMatrix.initH();
                    RealMatrix A = DegreeOfConditioningSymmetricMatrix.initA(L, H);

                    double Lmax = DegreeOfConditioningSymmetricMatrix.getMaxL(L);//известное
                    double Lmin = DegreeOfConditioningSymmetricMatrix.getMinL(L, eps);//известное
                    double m1 = DegreeOfConditioningSymmetricMatrix.conditionNumberOfTheMatrix(Lmin, Lmax);//известное число обусловленности

                    double Ln = DegreeOfConditioningSymmetricMatrix.directIterationMethod(A);//максимальное по модулю собственное значение
                    double L1 = DegreeOfConditioningSymmetricMatrix.reverseIterationMethod(A);//минимальное по модулю собственное значение
                    double m2 = DegreeOfConditioningSymmetricMatrix.conditionNumberOfTheMatrix(L1, Ln);//число обусловленности

                    /*absL1 = absL1 + Math.abs(Lmin - L1) / 10;
                    absLn = absLn + Math.abs(Lmax - Ln) / 10;
                    absm = absm + Math.abs(m1 - m2) / 10;*/

                    absL1 = Lmin - L1;//абсолютная погрешность = точное - приближенное
                    absLn =Lmax - Ln;
                    absM = m1 - m2;

                    otnosL1 = otnosL1 + Math.abs(absL1 / Lmin) * 100; //относительная погрешность в % = |абсолют / точное| * 100%
                    otnosLn = otnosLn + Math.abs(absLn / Lmax) * 100;
                    otnosM = otnosM + Math.abs(absM / m1) * 100;
                }
               /* System.out.println("N = " + n + ", L from " +  -eps + " to " + eps + ", " +
                        " accuracy L1 = " + Precision.round(absL1, 3) +
                        " accuracy Ln = " + Precision.round(absLn, 3) +
                        " accuracy m = " + Precision.round(absM, 3));*/

                System.out.println("N = " + n + ", L from " +  -eps + " to " + eps + ", " +
                        " средняя относительная погрешность в % для L1 = " + otnosL1 / 10+
                        " для Ln = " + otnosLn / 10+
                        " для m = " + otnosM / 10);
            }
        }
        System.exit(IER);
    }
}
