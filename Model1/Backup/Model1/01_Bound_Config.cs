using System;

public partial class globalmembers
{
    public const int SRC1 = 01;
    public const int SNK_0BAR = 02;



    public static void Bound_config()

    {
        // SRC1
        boundary[SRC1].p = 5; // BAR
        boundary[SRC1].t = 90; // DEGREE C
        boundary[SRC1].outlet[0] = 1; // OUTLET STREAM NO.
        boundary[SRC1].calc(0, 1);
        //
        // SNK_0BAR
        boundary[SNK_0BAR].p = 0.2F; // BAR
        boundary[SNK_0BAR].inlet[0] = 2;
        boundary[SNK_0BAR].inlet[1] = 3;
        boundary[SNK_0BAR].calc(2, 0);
        //   

    
    }

}
