using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Timers;

public partial class globalmembers
{
    internal static class Arrays
    {
        public static T[] InitializeWithDefaultInstances<T>(int length) where T : new()
        {
            T[] array = new T[length];
            for (int i = 0; i < length; i++)
            {
                array[i] = new T();
            }
            return array;
        }

        public static void DeleteArray<T>(T[] array) where T : System.IDisposable
        {
            foreach (T element in array)
            {
                if (element != null)
                    element.Dispose();
            }
        }
    }



    public class global_variable
    {
        public int sim_spd;
        public int run_count;
    }


    public class VALVE
    {
        // config
        public float cv;
        public float timeop; // time to open 100% in sec
        public float timecl; // time to close 100% in sec

        public string tagname = new string(new char[25]);
        public string description = new string(new char[50]);

        // calculated
        public float op; // command
        public double @out; // output
    }

    public class flow
    {

        public double massf; // kg/s
        public double p_in; // kg/m2
        public double p_out;
        public double t_in; //C
        public double t_out;
        //float density; // kg/m3

        public double k; // calculated flow conductance
        public double cp;
        public float h;
        public int rev_status;
        //////////////////////////////////////

        public double molwt;
        public double molf; // Kmol/hr
        public double ro; // mixed phase density Kg/m3

    }



    public class mechstream
    {

        public double pow; // kw
        public float speed; // rpm

    }








    public class HDR
    {
        public int[] inlet = new int[30];
        public int[] outlet = new int[30];
        public double p;
        public double t;
        public float vol;
        public float mass_in;
        public float mass_out;
        public float density;
        public float cp;
        public float ua;
        public float mass_content_of_hdr;
        //////////////////////////

        public float molwt;
        public double mvol;


        public void calc(int @in, int @out)
        {
            //calc_pressure(in,out);
            calc_c(@in, @out);
            calc_p(@in, @out);
            calc_temp(@in, @out);


            /////////////// write properties to stream ///////////////////////////////
            for (int x = 0; x < @out; x++)
            {
                stream[outlet[x]].p_in = p;
                if (stream[outlet[x]].rev_status == 1)
                {
                    stream[outlet[x]].t_in = t;
                    stream[outlet[x]].ro = density;
                    stream[outlet[x]].cp = cp;
                    stream[outlet[x]].molwt = molwt;

                }
            }
            for (int x = 0; x < @in; x++)
            {
                stream[inlet[x]].p_out = p;
                if (stream[inlet[x]].rev_status == -1)
                {
                    stream[inlet[x]].t_out = t;
                    stream[inlet[x]].ro = density;
                    stream[inlet[x]].cp = cp;
                    stream[outlet[x]].molwt = molwt;
                }
            }


        }

        public void calc_temp(int @in, int @out)
        {


            //////CALCULATE DENSITY OF MIXTURE IN HEADER////////////

            density = 1.0F; // g/cm3


            /////////////////////////calculate cp and temperature/////////////////////////////////////////////////////

            double sum_m_t;
            double sum_m;

            sum_m_t = 0F;
            sum_m = 0F;

            for (int x = 0; x < @in; x++)
            {
                if (stream[inlet[x]].rev_status == 1)
                {
                    sum_m_t += stream[inlet[x]].massf * stream[inlet[x]].t_out;
                    sum_m += stream[inlet[x]].massf;
                }
            }
            for (int x = 0; x < @out; x++)
            {
                if (stream[outlet[x]].rev_status == -1)
                {
                    sum_m_t += -stream[outlet[x]].massf * stream[outlet[x]].t_in;
                    sum_m += -stream[outlet[x]].massf;
                }
            }
            cp = 1.0F; //  (sum_m_cp + (mass_content_of_hdr *cp))/(mass_content_of_hdr+sum_mass+0.000000001); // cp calculated

            sum_m_t = sum_m_t + 10 * t;
            sum_m = sum_m + 10;

            // temperature caalculate /// need to add heat loss to ambient and the n final temperature
            t += (sum_m_t / sum_m) - t;

            // t= (sum_m_cp_t / (sum_m_cp+0.00000001));

            ////heat loss to ambient ////

            //      float heat_loss;

            //ua = 0.01; // j/m .s. K

            //heat_loss = ua * (t - TAMBIENT)* TIME_STEP ;

            // t = (sum_m_cp_t - heat_loss)/(sum_m_cp+0.00000001) ;

        }

        public void calc_c(int @in, int @out)
        {

            //  float no_of_component;
            float[] moles = new float[20];
            double sum_moles_in;
            double sum_moles_out;
            //   float sum_mf;





            sum_moles_in = 0F;
            for (int x = 0; x < @in; x++)
            {
                if (stream[inlet[x]].rev_status == 1)
                {
                    sum_moles_in += stream[inlet[x]].molf;
                }
            }
            mvol += sum_moles_in;

            molwt = 18.01528F;


            sum_moles_out = 0F;
            for (int x = 0; x < @out; x++)
            {
                if (stream[outlet[x]].rev_status == 1)
                {
                    sum_moles_out += stream[outlet[x]].molf;
                }
            }
            mvol -= sum_moles_out;


        }

        public void calc_p(int @in, int @out)
        {
            //8.31447;  // Gas Constant, R = 8.3145	kPa.m3 / ( Kg-mol.K)	
            p = (mvol * 8.31447 * (273.15 + t) / vol) / 100000; // pa to bar R = 8.31447

        }
        public HDR()
        {
            p = 0;
            //output=0;
        }

    }

    public class BOUNDARY
    {
        public int[] inlet = new int[30];
        public int[] outlet = new int[30];
        public float p; // bar
        public float t; // degree C
        public float molwt; // molecular wt
        public double ro; // mixed phase density Kg/m3
        public double h; // enthalpy kj/kg
        public double vf;



        public void calc(int @in, int @out)
        {

            h = st1.h_pT(p, t); // kj/kg
            ro = st1.rho_pT(p, t); // Kg/m3
            vf = st1.x_ph(p, h);
            double cp;
            cp = st1.Cp_pT(p, t); // kJ/(kg°C)


            ///////////////// CALCULATE MOL WT
            molwt = 18.01528F; // kg/kg-mol
            ///////////////////////////////////


            ///////////////passing values to outlet streams//////////////////////////////////////////////////
            for (int x = 0; x < @in; x++)
            {
                stream[inlet[x]].p_out = p;
            }

            for (int x = 0; x < @out; x++)
            {
                stream[outlet[x]].p_in = p;
                stream[outlet[x]].t_in = t;
                stream[outlet[x]].cp = cp;
                stream[outlet[x]].molwt = molwt;
                stream[outlet[x]].ro = ro;

            }
        }


        public BOUNDARY()
        {
            //
            p = 0;
        }


    }

    public class CONTROLLER
    {

        public float sp;
        public float mv;
        public float pv;
        public double vp;
        public double @out;
        public byte am;
        public byte lr;
        public byte rd;
        public byte @lock;

        public float p;
        public float i;
        public float d;
        public double iop;
        public float dop;

        public float mvr1;
        public float mvr2;
        public float err;


        public void calc()
        {
            float ee;
            float error;

            if (mvr1 == mvr2)
            {
                return;
            }
            if (rd == 1) // direct 1    // reverse 0
            {
                error = (mv - sp) / (mvr2 - mvr1) * 100;
            }
            else
            {
                error = (sp - mv) / (mvr2 - mvr1) * 100;
            }

            ee = err; // Assigning Previous cycle error

            if (am == 1)
            {
                iop = @out;
                vp = @out;
                err = error;
                return;
            }
            if (i <= (float)0.0)
            {
                i = (float)0.01;
            }

            iop = iop + p * error * 1 * 0.0166 / i;

            if (iop < 0F)
            {
                iop = 0F;
            }
            if (iop > 100F)
            {
                iop = 100F;
            }

            dop = p * d * 60 * (error - ee) / 1;


            if (@lock == 0)
            {
                @out = (p * error) + iop + dop;
            }

            if (@out < 0.01F)
            {
                @out = 0F;
            }
            if (@out > 99.99F)
            {
                @out = 100F;
            }
            err = error;
            vp = @out;


        }
        public CONTROLLER()
        {
            @out = 0F;
        }


    }

    public class HXSIMPLE
    {

        public float s_k; // shell and tube CV
        public float t_k;
        public float u;
        public float area;



        public void calc(int s_in, int t_in) // ( shell side stream inlet index no, tube side stream inlet index no)
        {

            double p_in;
            double p_out;
            double massf;
            float rho;
            int reverse;

            // calculate tube side flow//////////////////////////////////
            stream[t_in].k = t_k;

            p_in = stream[t_in].p_in;
            p_out = stream[t_in].p_out;
            massf = stream[t_in].massf;

            rho = 1F; // stream[indexno].density/1000 ;
            reverse = 1;
            if (p_in < p_out)
            {
                reverse = -1;
            }
            massf = ((t_k * Math.Sqrt(Math.Abs(p_in - p_out)) * reverse * rho));

            stream[t_in].massf = massf;
            stream[t_in].rev_status = reverse;
            // calculate shell side flow//////////////////////////////////
            stream[s_in].k = s_k;

            p_in = stream[s_in].p_in;
            p_out = stream[s_in].p_out;
            massf = stream[s_in].massf;

            rho = 1F; // stream[indexno].density/1000 ;
            reverse = 1;
            if (p_in < p_out)
            {
                reverse = -1;
            }
            massf = ((s_k * Math.Sqrt(Math.Abs(p_in - p_out)) * reverse * rho));

            stream[s_in].massf = massf;
            stream[s_in].rev_status = reverse;

            ////// calculate outlet temperature of shell nad tube

            double f_shell;
            double f_tube;
            double t_sin;
            double t_tin;
            double t_sout;
            double t_tout;
            double cp_t;
            double cp_s;
            double dt1;
            double dt2;
            double lmtd;
            double heat_acc_tube;
            double heat_acc_shell;
            double ratio;

            f_shell = stream[s_in].massf;
            f_tube = stream[t_in].massf;
            t_sin = stream[s_in].t_in;
            t_tin = stream[t_in].t_in;
            t_sout = stream[s_in].t_out;
            t_tout = stream[t_in].t_out;
            cp_s = stream[s_in].cp;
            cp_t = stream[t_in].cp;


            if (f_tube == 0F)
            {
                t_sout += 0.1 * (t_sin - t_sout);
                t_tout -= 0.001 * (t_tout - TAMBIENT);
            }

            if (f_shell == 0F)
            {
                t_tout -= 0.1 * (t_tout - t_tin);
                t_sout -= 0.001 * (t_sout - TAMBIENT);
            }

            if (f_tube == 0F && f_shell == 0F)
            {
                t_tout -= 0.001 * (t_tout - TAMBIENT);
                t_sout -= 0.001 * (t_sout - TAMBIENT);
            }

            if (f_tube != 0F && f_shell != 0F)
            {
                dt1 = t_sout - t_tin;
                dt2 = t_sin - t_tout;
                ratio = dt1 / dt2;
                if (ratio < 0.001F)
                {
                    lmtd = 0F; //(dt1+dt2)/2.0;
                }
                else
                {
                    lmtd = (dt1 - dt2) / Math.Log(ratio);
                }
                if (lmtd <= 0F)
                {
                    lmtd = 0F;
                }
                if (lmtd >= 100F)
                {
                    lmtd = 100F;
                }

                heat_acc_tube = u * area * lmtd - f_tube * cp_t * (t_tout - t_tin);
                t_tout += heat_acc_tube / 3600;
                // range (32, &t_tout, 200);

                heat_acc_shell = -u * area * lmtd + f_shell * cp_s * (t_sin - t_sout);
                t_sout += heat_acc_shell / 3600;
                // range (32, &t_sout, 200);
            }

            stream[s_in].t_out = t_sout;
            stream[t_in].t_out = t_tout;


        }

        public HXSIMPLE()
        {
            return;

        }

    }

    public class PUMP
    {
        public int StreamNo;
        public int MStreamNum;
        public int sucValveNum;
        public int disValveNum;
        public float Qscale; // scale factor for Flow. m3/hr
        public float DHscale; // scale factor for Head. m
        public float ETAscale; // scale factor for efficiency. fraction
        public float SpeedRef; // reference speed
        public float Speed;
        public float J; //wind milling flow conductance
        public float Kjr; // reverse flow factor

        // TUNING PARAMETRES
        // hsexp // exponent head for speed
        // qsexp ;// exponent for q

        // calculated values

        public double DH; // Actual Head m
        public float Eff; // PumpEfficiency fraction
        public double Pow; // pump Power KW
        public double Q; // Actual vol flow m3/sec



        public void calc()
        {

            //The available Head is calculated based on the differential pressure across the Pump. 
            // DH = 1000*deltaP/9.81*Rf*MW
            //
            //Rf- inlet stream density kg-mol/m3
            //deltaP - KPa
            //MW- Kg/kg-mol

            double dp = ((stream[StreamNo].p_out - stream[StreamNo].p_in)) * 100; // *100 bar to KPa
            float Rf = 1000.0F; // float Rf = stream[StreamNo].ro; // kg/m3

            DH = 1000 * dp / (9.81 * Rf);

            // default Pump Curve Q = f(DH) //  q = -1.3018 dh^3 + 1.7823 dh^2 -1.3108 dh + 1.8308
            Speed = mstream[MStreamNum].speed;
            float sp = Speed / SpeedRef;
            double dhn = DH / (DHscale * (sp * sp) + 0.00000001);

            double qn = (-1.3018 * dhn * dhn * dhn) + (1.7823 * dhn * dhn) - (1.3108 * dhn) + 1.8308;

            Q = qn * Qscale * sp; // * 3600.0; // m3/hr
            Q = Q * valve[sucValveNum].@out * valve[disValveNum].@out;
            if (Q <= 0F)
            {
                Q = 0.0F;
            }

            double massf = Q * Rf;

            double molf = massf / stream[StreamNo].molwt;
            //////////////////////////////////////////////////////////////////////////////////

            /// pump power required

            double Hpow = (Q / 3600.0) * DH * Rf * 9.81 / 1000.0; // hydraulic power KW

            /////////////////////////// calculate pump efficiency from curve. Eff right now constant efficiency 75%

            Eff = 0.75F;
            Pow = Hpow / Eff; // pump shaft Power

            /////////////////////////////////////////////////////////////////////////////////////
            //change in enthalpy, Tempchange  

            double deltaH = Pow / (molf + 0.0000001);
            if (Q <= 0.05F)
            {
                deltaH = 0.0F;
            }
            double t_out = (deltaH / stream[StreamNo].cp) + stream[StreamNo].t_in;

            ////////////////////////////////////////////////////////////////////////////////////

            stream[StreamNo].rev_status = 1; // means no reverxe flow;
            stream[StreamNo].t_out = t_out;
            stream[StreamNo].massf = massf;
            stream[StreamNo].molf = molf;
            mstream[MStreamNum].pow = Pow;
        }



        public PUMP()
        {
            //
            Q = 0F;

        }


    }

    public class SIMPLEDRUM
    {
        public int[] inlet = new int[30];
        public int[] v_out = new int[30];
        public int[] l_out = new int[30];

        public float elevation; // m
        public float orientation; // horizontal = 0 , vertical = 1
        public float dia; // diameter m
        public float height; // m
        public double mvol_liq; // lig mols // initialise
        public double mvol_vap; // vapour ,ols // initialise



        // calculated
        public double p;
        public double t;
        public double lvl;



        public void calc(int @in, int lout, int vout)
        {

            double vol;
            double sum_liq_mols_in;
            double sum_liq_mols_out;
            double vf;
            double liqVol;
            double vapVol;

            vol = (PI * dia * dia / 4) * height; //m3


            sum_liq_mols_in = 0F;
            for (int x = 0; x < @in; x++)
            {
                if (stream[inlet[x]].rev_status == 1)
                {
                    vf = 0.01F; // calculate vf for stream vf = stream[inlet[x]].vf ;
                    sum_liq_mols_in += stream[inlet[x]].molf * (1 - vf) / 3600; // kg/hr to kg/sec
                }
            }
            mvol_liq += sum_liq_mols_in;

            sum_liq_mols_out = 0F;
            for (int x = 0; x < lout; x++)
            {
                if (stream[l_out[x]].rev_status == 1)
                {
                    sum_liq_mols_out += stream[l_out[x]].molf / 3600; // kg/hr to kg/sec
                }
            }
            mvol_liq -= sum_liq_mols_out;

            // liqVol = mvol_liq* R* (t + 273.15)/(p * 100000) ;   // p =(mvol*R*(273.15+t)/vol)/100000; 
            liqVol = mvol_liq * R * (t + 273.15) / (1 * 100000);

            lvl = liqVol / vol * height;
            if (lvl > height)
            {
                lvl = height;
            }


            vapVol = (PI * dia * dia / 4) * ((1.1 * height) - lvl);

            double sum_vap_mols_in;
            double sum_vap_mols_out;

            sum_vap_mols_in = 0F;
            for (int x = 0; x < @in; x++)
            {
                if (stream[inlet[x]].rev_status == 1)
                {
                    vf = 0.01F; // calculate vf for stream vf = stream[inlet[x]].vf ;
                    sum_vap_mols_in += stream[inlet[x]].molf * vf / 3600;
                }
            }
            mvol_vap += sum_vap_mols_in;

            sum_vap_mols_out = 0F;
            for (int x = 0; x < vout; x++)
            {
                if (stream[v_out[x]].rev_status == 1)
                {
                    sum_vap_mols_out += stream[v_out[x]].molf / 3600;
                }
            }
            mvol_vap -= sum_vap_mols_out;


            // vapour pressure // drum pressure 

            p = mvol_vap * R * (t + 273.15) / (vapVol * 100000); // p bar // // p =(mvol*R*(273.15+t)/vol)/100000;

            // preesure at liquid bottol 
            double p_bottom;

            p_bottom = p + ((lvl + elevation) * 9.81 * 1000) / 100000; // ro = 1000 kg/m3  p in bar


            // pressure at Inlet stream
            // float p_inlet_stream ;
            //  p_inlet_stream = p + ((max ((lvl - port_height), 0 ) + elevation )  *9.81* 1000)/100000   ; 

            // calculate tmperature of drum 



            double sum_m_t;
            double sum_m;

            sum_m_t = 0F;
            sum_m = 0F;

            for (int x = 0; x < @in; x++)
            {
                if (stream[inlet[x]].rev_status == 1)
                {
                    sum_m_t += stream[inlet[x]].massf * stream[inlet[x]].t_out;
                    sum_m += stream[inlet[x]].massf;
                }
            }

            sum_m_t = sum_m_t + 10 * t; //  10 = ( mvol_liq + mvol_vap )* mol wt
            sum_m = sum_m + 10;

            // temperature caalculate /// need to add heat loss to ambient and the n final temperature
            t += (sum_m_t / sum_m) - t;

            // t= (sum_m_cp_t / (sum_m_cp+0.00000001));

            ////heat loss to ambient ////

    //        float heat_loss;

            //ua = 0.01; // j/m .s. K

            //heat_loss = ua * (t - TAMBIENT)* TIME_STEP ;

            // t = (sum_m_cp_t - heat_loss)/(sum_m_cp+0.00000001) ;




            ///////////////////////////////////////////////////////////////////////////////////////////////

            /////////////// write properties to stream ///////////////////////////////
            for (int x = 0; x < lout; x++)
            {
                stream[l_out[x]].p_in = p_bottom;
                if (stream[l_out[x]].rev_status == 1)
                {
                    stream[l_out[x]].t_in = t;
                    stream[l_out[x]].ro = 1; // density =1 gm/cm3
                    stream[l_out[x]].cp = 1; // cp = 1
                    stream[l_out[x]].molwt = 18.01528;
                }
            }
            //////////////////////////////////////////////
            for (int x = 0; x < vout; x++)
            {
                stream[v_out[x]].p_in = p;
                if (stream[v_out[x]].rev_status == 1)
                {
                    stream[v_out[x]].t_in = t;
                    stream[v_out[x]].ro = 1; // density =1 gm/cm3
                    stream[v_out[x]].cp = 1; // cp = 1
                    stream[v_out[x]].molwt = 18.01528;
                }
            }
            //////////////////////////////////////////////

            for (int x = 0; x < @in; x++)
            {
                stream[inlet[x]].p_out = p; //     //  p_inlet_stream
                if (stream[inlet[x]].rev_status == -1)
                {
                    stream[inlet[x]].t_out = t;
                    stream[inlet[x]].ro = 1;
                    stream[inlet[x]].cp = 1;
                    stream[inlet[x]].molwt = 18.01528;
                }
            }








        }



        public SIMPLEDRUM()
        {
            //
            p = 0;

        }


    }

    ////////////////////////////////////////////////////////////////////////////functions from cpp////////////////////////////////////////////////

    public void calc_molar_flow()
    {
        int i;
        for (i = 0; i < 1000; i++)
        {
            stream[i].molf = stream[i].massf / stream[i].molwt;
        }
    }



    public static void calc_valve(int indexno)
    {
        float op;
        double @out;
        float timeop;
        float timecl;

        op = valve[indexno].op;
        @out = valve[indexno].@out;
        timeop = valve[indexno].timeop;
        timecl = valve[indexno].timecl;

        if (op > @out)
        {
            @out += TIME_STEP / timeop;
            if (@out > op)
            {
                @out = op;
            }
        }
        if (op < @out)
        {
            @out += -TIME_STEP / timecl;
            if (@out < op)
            {
                @out = op;
            }
        }

        valve[indexno].@out = @out;

    }

    public static void calc_stream(int indexno, int valve_indexno)
    {
        calc_valve(valve_indexno);

        double p_in;
        double p_out;
        double t_in;
        double t_out;
        double massf;
     //   float h;
        double cv;
        double ro; //kg/m3
        int reverse;
        float trim;
        trim = 1F;
        double j; // flow conductance

        cv = valve[valve_indexno].cv * valve[valve_indexno].@out * trim;

        // Flow conductance can be calculated from C v  with the following equation. j = 	0.00075379*cv

        j = 0.00075379 * cv; //( kg/hr / sqrt(Pa.Kg/m3)  )


        stream[indexno].k = j;



        p_in = stream[indexno].p_in * 100000; // convert to pa
        p_out = stream[indexno].p_out * 100000; // convert to pa
        massf = stream[indexno].massf; // kg/hr
        t_in = stream[indexno].t_in;
        t_out = stream[indexno].t_out;
        ro = stream[indexno].ro; // kg/m3

        reverse = 1;
        if (p_in < p_out)
        {
            reverse = -1;
        }


        //
        massf = reverse * j * Math.Sqrt(Math.Abs(p_in - p_out) * ro); //kg/hr


        //	massf += ((cv * sqrt(abs(p_in-p_out)) *reverse* rho)-massf)*TIME_STEP;

        if (reverse == 1 )
        {
            t_out = t_in;
        }
        else
        {
            t_in = t_out;
        }


        // calculate Cp 
        // cp_calc

        //	stream[indexno].cp = 1 ;
        //cp = stream[indexno].cp ;

        // calculate enthalpy 

        // 	  h = massf *  stream[indexno].cp * t_in * reverse;


        stream[indexno].massf = massf;
        stream[indexno].rev_status = reverse;
        stream[indexno].t_in = t_in;
        stream[indexno].t_out = t_out;

    }



}
