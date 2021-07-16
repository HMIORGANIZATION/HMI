public partial class globalmembers
{
    public const int NO_OF_VLV = 100;
    //
    public const int FV001 = 01;
    public const int MAN_VLV_01 = 02;
    public const int MAN_VLV_02 = 03;

    public static void valve_config()
    {
        //FV001
        valve[FV001].cv = 20;
        valve[FV001].timeop = 5;
        valve[FV001].timecl = 5;
        valve[1].op = 0;
        //
        //MAN_VLV_01
        valve[MAN_VLV_01].cv = 10;
        valve[MAN_VLV_01].timeop = 1;
        valve[MAN_VLV_01].timecl = 1;
        //
        //MAN_VLV_02
        valve[MAN_VLV_02].cv = 10;
        valve[MAN_VLV_02].timeop = 1;
        valve[MAN_VLV_02].timecl = 1;
        //
    }

}