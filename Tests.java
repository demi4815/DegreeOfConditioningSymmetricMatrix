package vsu.labs;

import org.apache.commons.math3.linear.RealMatrix;

public class Tests
{
    public static void testInitA()
    {
        RealMatrix L = DegreeOfConditioningSymmetricMatrix.initL(1);
        RealMatrix H = DegreeOfConditioningSymmetricMatrix.initH();
        RealMatrix A = DegreeOfConditioningSymmetricMatrix.initA(L, H);

        System.out.println(L);
        System.out.println(H);
        System.out.println(A);
    }

    public static void testDirect()
    {
        RealMatrix L = DegreeOfConditioningSymmetricMatrix.initL(1);
        RealMatrix H = DegreeOfConditioningSymmetricMatrix.initH();
        RealMatrix A = DegreeOfConditioningSymmetricMatrix.initA(L, H);

        System.out.println(L);
        System.out.println(H);
        System.out.println(A);

        System.out.println();

        DegreeOfConditioningSymmetricMatrix.directIterationMethod(A);
        System.out.println();
        DegreeOfConditioningSymmetricMatrix.reverseIterationMethod(A);
    }
}
