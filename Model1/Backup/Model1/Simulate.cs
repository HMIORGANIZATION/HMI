
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Timers;





public partial class globalmembers
{
public const double TIME_STEP = 0.25;
public const double TAMBIENT = 35.0;
public const double PI = 3.14159265 ;
public const double R = 8.31447;  // Gas Constant, R = 8.3145	kPa.m3 / ( Kg-mol.K)	
//
public static steamtable st1 = new steamtable();
public static global_variable gl = new global_variable();
public static flow[] stream = Arrays.InitializeWithDefaultInstances<flow>(NO_OF_STREAMS);
public static VALVE[] valve = Arrays.InitializeWithDefaultInstances<VALVE>(NO_OF_VLV);
public static mechstream[] mstream = Arrays.InitializeWithDefaultInstances<mechstream>(1000);
public static HDR[] hdr = Arrays.InitializeWithDefaultInstances<HDR>(NO_OF_HDR);
public static BOUNDARY[] boundary = Arrays.InitializeWithDefaultInstances<BOUNDARY>(50);
public static CONTROLLER[] pid = Arrays.InitializeWithDefaultInstances<CONTROLLER>(500);
public static HXSIMPLE[] hxs = Arrays.InitializeWithDefaultInstances<HXSIMPLE>(500);
public static PUMP[] pump = Arrays.InitializeWithDefaultInstances<PUMP>(100);
public static SIMPLEDRUM[] sdrum = Arrays.InitializeWithDefaultInstances<SIMPLEDRUM>(200);
  
         

    public static void Simulation()
    {
        gl.run_count++;

        //// TEST MODEL ////////////////////////////////////////////////////////
        calc_stream(1, FV001);
        calc_stream(2,MAN_VLV_01);
        calc_stream(3, MAN_VLV_01);
        hdr[SUCTION_HDR].calc(1, 2);// INLET STREAM 1 ; OUTLET STREAMS 2 , 3 //
        ////////////////////////////////////////////////////////////////////////
        
        
      
    }
}
// }

