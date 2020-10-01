using System;
using System.Diagnostics;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Timers;





namespace Model1
{
    public class program

    {
   
             public static void Main()
            {
                var watch = new System.Diagnostics.Stopwatch();
                float TIMESTEP = 250; // 250 MILISECONDS = 0.25 SEC
                float simtime, elaps_time; 
                 simtime = 0;
                 elaps_time = 0;
                 bool isRunning = true;
                  globalmembers.gl.sim_spd = 100;

            globalmembers.model_config();


            while (isRunning)	
                 {

                        watch.Start();
                        while (elaps_time < TIMESTEP)
                        {
                                   if (globalmembers.gl.run_count < globalmembers.gl.sim_spd)
                                            globalmembers.Simulation();
                                   elaps_time = watch.ElapsedMilliseconds;
                         
                        }
                        watch.Stop();
                        simtime += (elaps_time * globalmembers.gl.run_count)/1000; // sec

                Console.Write("\r");
                //Console.Write(simtime);
                ///   Console.Write(globalmembers.hdr[1].p);



                globalmembers.gl.run_count = 0;
                        watch.Reset();
                        elaps_time = 0;
                 }





           
             }

           
        
    }
}


