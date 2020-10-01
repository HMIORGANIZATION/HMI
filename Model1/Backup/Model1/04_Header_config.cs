public partial class globalmembers
{

    public const int NO_OF_HDR = 100;
    public const int SUCTION_HDR = 01;

    public static void hdr_config()
    {
        // initialise all hdrs

        for (int i = 0; i <= NO_OF_HDR; i++)
        {
            hdr[i].molwt = 18.01528F;
            hdr[i].p = 1;
            hdr[i].t = 50;
            hdr[i].density = 1000; // initialise
            hdr[i].mvol = 50; //Kmol
            hdr[1].cp = 4.2F;// initialise

        }

        /////////////////////////////////////////////////

        // SUCTION_HDR
        hdr[SUCTION_HDR].vol = 1;  // M3
        hdr[SUCTION_HDR].ua = 0;// 0.1; // HT COEFFIENT TO ATMOSPHERE
        hdr[SUCTION_HDR].inlet[0] = 1; /// STREAM #1 IS HEADER INLET
        hdr[SUCTION_HDR].outlet[0] = 2; ///
        hdr[SUCTION_HDR].outlet[1] = 3; ///
        ///////////////////////////////////////////////////////




    }

}