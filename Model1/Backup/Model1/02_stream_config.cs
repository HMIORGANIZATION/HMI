public partial class globalmembers
{

    public const int NO_OF_STREAMS = 100;

    public static void stream_config()
    {
        // auto initialise all streams for composition

        for (int i = 0; i <= NO_OF_STREAMS; i++) // no of stream + 1 should be added here
        {

            stream[i].molwt = 18.01528;
            stream[i].cp = 4.2;
            stream[i].t_in = TAMBIENT;
            stream[i].t_out = TAMBIENT;
            stream[i].massf = 0.05;
            stream[i].molf = 0.05;
            stream[i].ro = 1000.0;
        }



    }

}