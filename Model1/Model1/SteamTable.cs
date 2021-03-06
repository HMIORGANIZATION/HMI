﻿using System;


public partial class globalmembers
{
    public class steamtable
    {


        //***********************************************************************************************************
        //*1.2 Temperature
        public double Tsat_p(double p)
        {
            double tempTsat_p = 0;
            p = toSIunit_p(p);
            if (p >= 0.000611657 && p <= 22.06395 + 0.001) //0.001 Added to enable the tripple point.
            {
                tempTsat_p = fromSIunit_T(T4_p(p));
            }
            else
            {
                tempTsat_p = 9999999999;  // error // need to check how to handle this
            }
            return tempTsat_p;
        }

        public double Tsat_s(double s)
        {
            double tempTsat_s = 0;
            s = toSIunit_s(s);
            if (s > -0.0001545495919 && s < 9.155759395)
            {
                tempTsat_s = fromSIunit_T(T4_p(p4_s(s)));
            }
            else
            {
                tempTsat_s = 9999999999;  // error // need to check how to handle this
            }
            return tempTsat_s;
        }

        public double T_ph(double p, double h)
        {
            double tempT_ph = 0;
            p = toSIunit_p(p);
            h = toSIunit_h(h);
            switch (region_ph(p, h))
            {
                case 1:
                    tempT_ph = fromSIunit_T(T1_ph(p, h));
                    break;
                case 2:
                    tempT_ph = fromSIunit_T(T2_ph(p, h));
                    break;
                case 3:
                    tempT_ph = fromSIunit_T(T3_ph(p, h));
                    break;
                case 4:
                    tempT_ph = fromSIunit_T(T4_p(p));
                    break;
                case 5:
                    tempT_ph = fromSIunit_T(T5_ph(p, h));
                    break;
                default:
                    tempT_ph = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempT_ph;
        }

        public double T_ps(double p, double s)
        {
            double tempT_ps = 0;
            p = toSIunit_p(p);
            s = toSIunit_s(s);
            switch (region_ps(p, s))
            {
                case 1:
                    tempT_ps = fromSIunit_T(T1_ps(p, s));
                    break;
                case 2:
                    tempT_ps = fromSIunit_T(T2_ps(p, s));
                    break;
                case 3:
                    tempT_ps = fromSIunit_T(T3_ps(p, s));
                    break;
                case 4:
                    tempT_ps = fromSIunit_T(T4_p(p));
                    break;
                case 5:
                    tempT_ps = fromSIunit_T(T5_ps(p, s));
                    break;
                default:
                    tempT_ps = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempT_ps;
        }

        public double T_hs(double h, double s)
        {
            double tempT_hs = 0;
            h = toSIunit_h(h);
            s = toSIunit_s(s);
            switch (Region_hs(h, s))
            {
                case 1:
                    tempT_hs = fromSIunit_T(T1_ph(p1_hs(h, s), h));
                    break;
                case 2:
                    tempT_hs = fromSIunit_T(T2_ph(p2_hs(h, s), h));
                    break;
                case 3:
                    tempT_hs = fromSIunit_T(T3_ph(p3_hs(h, s), h));
                    break;
                case 4:
                    tempT_hs = fromSIunit_T(T4_hs(h, s));
                    break;
                case 5:
                    tempT_hs = 9999999999;  // error // need to check how to handle this //Functions of hs is not implemented in region 5
                    break;
                default:
                    tempT_hs = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempT_hs;
        }


        //***********************************************************************************************************
        //*1.3 Pressure (p)
        public double psat_T(double T)
        {
            double temppsat_T = 0;
            T = toSIunit_T(T);
            if (T <= 647.096 && T > 273.15)
            {
                temppsat_T = fromSIunit_p(p4_T(T));
            }
            else
            {
                temppsat_T = 9999999999;  // error // need to check how to handle this
            }
            return temppsat_T;
        }

        public double psat_s(double s)
        {
            double temppsat_s = 0;
            s = toSIunit_s(s);
            if (s > -0.0001545495919 && s < 9.155759395)
            {
                temppsat_s = fromSIunit_p(p4_s(s));
            }
            else
            {
                temppsat_s = 9999999999;  // error // need to check how to handle this
            }
            return temppsat_s;
        }

        public double p_hs(double h, double s)
        {
            double tempp_hs = 0;
            h = toSIunit_h(h);
            s = toSIunit_s(s);
            switch (Region_hs(h, s))
            {
                case 1:
                    tempp_hs = fromSIunit_p(p1_hs(h, s));
                    break;
                case 2:
                    tempp_hs = fromSIunit_p(p2_hs(h, s));
                    break;
                case 3:
                    tempp_hs = fromSIunit_p(p3_hs(h, s));
                    break;
                case 4:
                    tempp_hs = fromSIunit_p(p4_T(T4_hs(h, s)));
                    break;
                case 5:
                    tempp_hs = 9999999999;  // error // need to check how to handle this //Functions of hs is not implemented in region 5
                    break;
                default:
                    tempp_hs = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempp_hs;
        }

        public double p_hrho(double h, double rho)
        {
            //*Not valid for water or sumpercritical since water rho does not change very much with p.
            //*Uses iteration to find p.
            double High_Bound = 0;
            double Low_Bound = 0;
            double p = 0;
            double rhos = 0;
            High_Bound = fromSIunit_p(100);
            Low_Bound = fromSIunit_p(0.000611657);
            p = fromSIunit_p(10);
            rhos = 1 / v_ph(p, h);
            while (Math.Abs(rho - rhos) > 0.0000001) // Math.Abs(Double)
            {
                rhos = 1 / v_ph(p, h);
                if (rhos >= rho)
                {
                    High_Bound = p;
                }
                else
                {
                    Low_Bound = p;
                }
                p = (Low_Bound + High_Bound) / 2;
            }
            return p;
        }
        //***********************************************************************************************************
        //*1.4 Enthalpy (h)

        public double hV_p(double p)
        {
            double temphV_p = 0;
            p = toSIunit_p(p);
            if (p > 0.000611657 && p < 22.06395)
            {
                temphV_p = fromSIunit_h(h4V_p(p));
            }
            else
            {
                temphV_p = 9999999999;  // error // need to check how to handle this
            }
            return temphV_p;
        }

        public double hL_p(double p)
        {
            double temphL_p = 0;
            p = toSIunit_p(p);
            if (p > 0.000611657 && p < 22.06395)
            {
                temphL_p = fromSIunit_h(h4L_p(p));
            }
            else
            {
                temphL_p = 9999999999;  // error // need to check how to handle this
            }
            return temphL_p;
        }

        public double hV_T(double T)
        {
            double temphV_T = 0;
            T = toSIunit_T(T);
            if (T > 273.15 && T < 647.096)
            {
                temphV_T = fromSIunit_h(h4V_p(p4_T(T)));
            }
            else
            {
                temphV_T = 9999999999;  // error // need to check how to handle this
            }
            return temphV_T;
        }

        public double hL_T(double T)
        {
            double temphL_T = 0;
            T = toSIunit_T(T);
            if (T > 273.15 && T < 647.096)
            {
                temphL_T = fromSIunit_h(h4L_p(p4_T(T)));
            }
            else
            {
                temphL_T = 9999999999;  // error // need to check how to handle this
            }
            return temphL_T;
        }

        public double h_pT(double p, double T)
        {
            double temph_pT = 0;
            p = toSIunit_p(p);
            T = toSIunit_T(T);
            switch (region_pT(p, T))
            {
                case 1:
                    temph_pT = fromSIunit_h(h1_pT(p, T));
                    break;
                case 2:
                    temph_pT = fromSIunit_h(h2_pT(p, T));
                    break;
                case 3:
                    temph_pT = fromSIunit_h(h3_pT(p, T));
                    break;
                case 4:
                    temph_pT = 9999999999;  // error // need to check how to handle this
                    break;
                case 5:
                    temph_pT = fromSIunit_h(h5_pT(p, T));
                    break;
                default:
                    temph_pT = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return temph_pT;
        }

        public double h_ps(double p, double s)
        {
            double temph_ps = 0;
            double xs = 0;
            p = toSIunit_p(p);
            s = toSIunit_s(s);
            switch (region_ps(p, s))
            {
                case 1:
                    temph_ps = fromSIunit_h(h1_pT(p, T1_ps(p, s)));
                    break;
                case 2:
                    temph_ps = fromSIunit_h(h2_pT(p, T2_ps(p, s)));
                    break;
                case 3:
                    temph_ps = fromSIunit_h(h3_rhoT(1 / v3_ps(p, s), T3_ps(p, s)));
                    break;
                case 4:
                    xs = x4_ps(p, s);
                    temph_ps = fromSIunit_h(xs * h4V_p(p) + (1 - xs) * h4L_p(p));
                    break;
                case 5:
                    temph_ps = fromSIunit_h(h5_pT(p, T5_ps(p, s)));
                    break;
                default:
                    temph_ps = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return temph_ps;
        }

        public double h_px(double p, double x)
        {
            double hL = 0;
            double hV = 0;
            p = toSIunit_p(p);
            x = toSIunit_x(x);
            if (x > 1 || x < 0 || p >= 22.064)
            {
                return 9999999999;  // error // need to check how to handle this
            }
            hL = h4L_p(p);
            hV = h4V_p(p);
            return fromSIunit_h(hL + x * (hV - hL));
        }
        public double h_Tx(double T, double x)
        {
            double hL = 0;
            double hV = 0;
            double p = 0;
            T = toSIunit_T(T);
            x = toSIunit_x(x);
            if (x > 1 || x < 0 || T >= 647.096)
            {
                return 9999999999;  // error // need to check how to handle this
            }
            p = p4_T(T);
            hL = h4L_p(p);
            hV = h4V_p(p);
            return fromSIunit_h(hL + x * (hV - hL));
        }

        public double h_prho(double p, double rho)
        {
            double temph_prho = 0;
            double hL = 0;
            double hV = 0;
            double vL = 0;
            double vV = 0;
            double x = 0;
            p = toSIunit_p(p);
            rho = 1 / toSIunit_v(1 / rho);
            switch (Region_prho(p, rho))
            {
                case 1:
                    temph_prho = fromSIunit_h(h1_pT(p, T1_prho(p, rho)));
                    break;
                case 2:
                    temph_prho = fromSIunit_h(h2_pT(p, T2_prho(p, rho)));
                    break;
                case 3:
                    temph_prho = fromSIunit_h(h3_rhoT(rho, T3_prho(p, rho)));
                    break;
                case 4:
                    if (p < 16.529)
                    {
                        vV = v2_pT(p, T4_p(p));
                        vL = v1_pT(p, T4_p(p));
                    }
                    else
                    {
                        vV = v3_ph(p, h4V_p(p));
                        vL = v3_ph(p, h4L_p(p));
                    }
                    hV = h4V_p(p);
                    hL = h4L_p(p);
                    x = (1 / rho - vL) / (vV - vL);
                    temph_prho = fromSIunit_h((1 - x) * hL + x * hV);
                    break;
                case 5:
                    temph_prho = fromSIunit_h(h5_pT(p, T5_prho(p, rho)));
                    break;
                default:
                    temph_prho = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return temph_prho;
        }

        //***********************************************************************************************************
        //*1.5 Specific Volume (v)
        public double vV_p(double p)
        {
            double tempvV_p = 0;
            p = toSIunit_p(p);
            if (p > 0.000611657 && p < 22.06395)
            {
                if (p < 16.529)
                {
                    tempvV_p = fromSIunit_v(v2_pT(p, T4_p(p)));
                }
                else
                {
                    tempvV_p = fromSIunit_v(v3_ph(p, h4V_p(p)));
                }
            }
            else
            {
                tempvV_p = 9999999999;  // error // need to check how to handle this
            }
            return tempvV_p;
        }

        public double vL_p(double p)
        {
            double tempvL_p = 0;
            p = toSIunit_p(p);
            if (p > 0.000611657 && p < 22.06395)
            {
                if (p < 16.529)
                {
                    tempvL_p = fromSIunit_v(v1_pT(p, T4_p(p)));
                }
                else
                {
                    tempvL_p = fromSIunit_v(v3_ph(p, h4L_p(p)));
                }
            }
            else
            {
                tempvL_p = 9999999999;  // error // need to check how to handle this
            }
            return tempvL_p;
        }

        public double vV_T(double T)
        {
            double tempvV_T = 0;
            T = toSIunit_T(T);
            if (T > 273.15 && T < 647.096)
            {
                if (T <= 623.15)
                {
                    tempvV_T = fromSIunit_v(v2_pT(p4_T(T), T));
                }
                else
                {
                    tempvV_T = fromSIunit_v(v3_ph(p4_T(T), h4V_p(p4_T(T))));
                }
            }
            else
            {
                tempvV_T = 9999999999;  // error // need to check how to handle this
            }
            return tempvV_T;
        }

        public double vL_T(double T)
        {
            double tempvL_T = 0;
            T = toSIunit_T(T);
            if (T > 273.15 && T < 647.096)
            {
                if (T <= 623.15)
                {
                    tempvL_T = fromSIunit_v(v1_pT(p4_T(T), T));
                }
                else
                {
                    tempvL_T = fromSIunit_v(v3_ph(p4_T(T), h4L_p(p4_T(T))));
                }
            }
            else
            {
                tempvL_T = 9999999999;  // error // need to check how to handle this
            }
            return tempvL_T;
        }

        public double v_pT(double p, double T)
        {
            double tempv_pT = 0;
            p = toSIunit_p(p);
            T = toSIunit_T(T);
            switch (region_pT(p, T))
            {
                case 1:
                    tempv_pT = fromSIunit_v(v1_pT(p, T));
                    break;
                case 2:
                    tempv_pT = fromSIunit_v(v2_pT(p, T));
                    break;
                case 3:
                    tempv_pT = fromSIunit_v(v3_ph(p, h3_pT(p, T)));
                    break;
                case 4:
                    tempv_pT = 9999999999;  // error // need to check how to handle this
                    break;
                case 5:
                    tempv_pT = fromSIunit_v(v5_pT(p, T));
                    break;
                default:
                    tempv_pT = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempv_pT;
        }
        public double v_ph(double p, double h)
        {
            double tempv_ph = 0;
            double xs = 0;
            double v4V = 0;
            double v4L = 0;
            p = toSIunit_p(p);
            h = toSIunit_h(h);
            switch (region_ph(p, h))
            {
                case 1:
                    tempv_ph = fromSIunit_v(v1_pT(p, T1_ph(p, h)));
                    break;
                case 2:
                    tempv_ph = fromSIunit_v(v2_pT(p, T2_ph(p, h)));
                    break;
                case 3:
                    tempv_ph = fromSIunit_v(v3_ph(p, h));
                    break;
                case 4:
                    xs = x4_ph(p, h);
                    if (p < 16.529)
                    {
                        v4V = v2_pT(p, T4_p(p));
                        v4L = v1_pT(p, T4_p(p));
                    }
                    else
                    {
                        v4V = v3_ph(p, h4V_p(p));
                        v4L = v3_ph(p, h4L_p(p));
                    }
                    tempv_ph = fromSIunit_v((xs * v4V + (1 - xs) * v4L));
                    break;
                case 5:
                    tempv_ph = fromSIunit_v(v5_pT(p, T5_ph(p, h)));
                    break;
                default:
                    tempv_ph = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempv_ph;
        }

        public double v_ps(double p, double s)
        {
            double tempv_ps = 0;
            double xs = 0;
            double v4V = 0;
            double v4L = 0;
            p = toSIunit_p(p);
            s = toSIunit_s(s);
            switch (region_ps(p, s))
            {
                case 1:
                    tempv_ps = fromSIunit_v(v1_pT(p, T1_ps(p, s)));
                    break;
                case 2:
                    tempv_ps = fromSIunit_v(v2_pT(p, T2_ps(p, s)));
                    break;
                case 3:
                    tempv_ps = fromSIunit_v(v3_ps(p, s));
                    break;
                case 4:
                    xs = x4_ps(p, s);
                    if (p < 16.529)
                    {
                        v4V = v2_pT(p, T4_p(p));
                        v4L = v1_pT(p, T4_p(p));
                    }
                    else
                    {
                        v4V = v3_ph(p, h4V_p(p));
                        v4L = v3_ph(p, h4L_p(p));
                    }
                    tempv_ps = fromSIunit_v((xs * v4V + (1 - xs) * v4L));
                    break;
                case 5:
                    tempv_ps = fromSIunit_v(v5_pT(p, T5_ps(p, s)));
                    break;
                default:
                    tempv_ps = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempv_ps;
        }
        //***********************************************************************************************************
        //*1.6 Density (rho)
        // Density is calculated as 1/v

        public double rhoV_p(double p)
        {
            return 1 / vV_p(p);
        }

        public double rhoL_p(double p)
        {
            return 1 / vL_p(p);
        }

        public double rhoL_T(double T)
        {
            return 1 / vL_T(T);
        }

        public double rhoV_T(double T)
        {
            return 1 / vV_T(T);
        }

        public double rho_pT(double p, double T)
        {
            return 1 / v_pT(p, T);
        }

        public double rho_ph(double p, double h)
        {
            return 1 / v_ph(p, h);
        }

        public double rho_ps(double p, double s)
        {
            return 1 / v_ps(p, s);
        }


        //***********************************************************************************************************
        //*1.7 Specific entropy (s)

        public double sV_p(double p)
        {
            double tempsV_p = 0;
            p = toSIunit_p(p);
            if (p > 0.000611657 && p < 22.06395)
            {
                if (p < 16.529)
                {
                    tempsV_p = fromSIunit_s(s2_pT(p, T4_p(p)));
                }
                else
                {
                    tempsV_p = fromSIunit_s(s3_rhoT(1 / (v3_ph(p, h4V_p(p))), T4_p(p)));
                }
            }
            else
            {
                tempsV_p = 9999999999;  // error // need to check how to handle this
            }
            return tempsV_p;
        }

        public double sL_p(double p)
        {
            double tempsL_p = 0;
            p = toSIunit_p(p);
            if (p > 0.000611657 && p < 22.06395)
            {
                if (p < 16.529)
                {
                    tempsL_p = fromSIunit_s(s1_pT(p, T4_p(p)));
                }
                else
                {
                    tempsL_p = fromSIunit_s(s3_rhoT(1 / (v3_ph(p, h4L_p(p))), T4_p(p)));
                }
            }
            else
            {
                tempsL_p = 9999999999;  // error // need to check how to handle this
            }
            return tempsL_p;
        }

        public double sV_T(double T)
        {
            double tempsV_T = 0;
            T = toSIunit_T(T);
            if (T > 273.15 && T < 647.096)
            {
                if (T <= 623.15)
                {
                    tempsV_T = fromSIunit_s(s2_pT(p4_T(T), T));
                }
                else
                {
                    tempsV_T = fromSIunit_s(s3_rhoT(1 / (v3_ph(p4_T(T), h4V_p(p4_T(T)))), T));
                }
            }
            else
            {
                tempsV_T = 9999999999;  // error // need to check how to handle this
            }
            return tempsV_T;
        }

        public double sL_T(double T)
        {
            double tempsL_T = 0;
            T = toSIunit_T(T);
            if (T > 273.15 && T < 647.096)
            {
                if (T <= 623.15)
                {
                    tempsL_T = fromSIunit_s(s1_pT(p4_T(T), T));
                }
                else
                {
                    tempsL_T = fromSIunit_s(s3_rhoT(1 / (v3_ph(p4_T(T), h4L_p(p4_T(T)))), T));
                }
            }
            else
            {
                tempsL_T = 9999999999;  // error // need to check how to handle this
            }
            return tempsL_T;
        }

        public double s_pT(double p, double T)
        {
            double temps_pT = 0;
            p = toSIunit_p(p);
            T = toSIunit_T(T);
            switch (region_pT(p, T))
            {
                case 1:
                    temps_pT = fromSIunit_s(s1_pT(p, T));
                    break;
                case 2:
                    temps_pT = fromSIunit_s(s2_pT(p, T));
                    break;
                case 3:
                    temps_pT = fromSIunit_s(s3_rhoT(1 / v3_ph(p, h3_pT(p, T)), T));
                    break;
                case 4:
                    temps_pT = 9999999999;  // error // need to check how to handle this
                    break;
                case 5:
                    temps_pT = fromSIunit_s(s5_pT(p, T));
                    break;
                default:
                    temps_pT = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return temps_pT;
        }

        public double s_ph(double p, double h)
        {
            double temps_ph = 0;
            double Ts = 0;
            double xs = 0;
            double s4V = 0;
            double s4L = 0;
            double v4V = 0;
            double v4L = 0;
            p = toSIunit_p(p);
            h = toSIunit_h(h);
            switch (region_ph(p, h))
            {
                case 1:
                    temps_ph = fromSIunit_s(s1_pT(p, T1_ph(p, h)));
                    break;
                case 2:
                    temps_ph = fromSIunit_s(s2_pT(p, T2_ph(p, h)));
                    break;
                case 3:
                    temps_ph = fromSIunit_s(s3_rhoT(1 / v3_ph(p, h), T3_ph(p, h)));
                    break;
                case 4:
                    Ts = T4_p(p);
                    xs = x4_ph(p, h);
                    if (p < 16.529)
                    {
                        s4V = s2_pT(p, Ts);
                        s4L = s1_pT(p, Ts);
                    }
                    else
                    {
                        v4V = v3_ph(p, h4V_p(p));
                        s4V = s3_rhoT(1 / v4V, Ts);
                        v4L = v3_ph(p, h4L_p(p));
                        s4L = s3_rhoT(1 / v4L, Ts);
                    }
                    temps_ph = fromSIunit_s((xs * s4V + (1 - xs) * s4L));
                    break;
                case 5:
                    temps_ph = fromSIunit_s(s5_pT(p, T5_ph(p, h)));
                    break;
                default:
                    temps_ph = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return temps_ph;
        }

        //***********************************************************************************************************
        //*1.8 Specific internal energy (u)
        public double uV_p(double p)
        {
            double tempuV_p = 0;
            p = toSIunit_p(p);
            if (p > 0.000611657 && p < 22.06395)
            {
                if (p < 16.529)
                {
                    tempuV_p = fromSIunit_u(u2_pT(p, T4_p(p)));
                }
                else
                {
                    tempuV_p = fromSIunit_u(u3_rhoT(1 / (v3_ph(p, h4V_p(p))), T4_p(p)));
                }
            }
            else
            {
                tempuV_p = 9999999999;  // error // need to check how to handle this
            }
            return tempuV_p;
        }

        public double uL_p(double p)
        {
            double tempuL_p = 0;
            p = toSIunit_p(p);
            if (p > 0.000611657 && p < 22.06395)
            {
                if (p < 16.529)
                {
                    tempuL_p = fromSIunit_u(u1_pT(p, T4_p(p)));
                }
                else
                {
                    tempuL_p = fromSIunit_u(u3_rhoT(1 / (v3_ph(p, h4L_p(p))), T4_p(p)));
                }
            }
            else
            {
                tempuL_p = 9999999999;  // error // need to check how to handle this
            }
            return tempuL_p;
        }

        public double uV_T(double T)
        {
            double tempuV_T = 0;
            T = toSIunit_T(T);
            if (T > 273.15 && T < 647.096)
            {
                if (T <= 623.15)
                {
                    tempuV_T = fromSIunit_u(u2_pT(p4_T(T), T));
                }
                else
                {
                    tempuV_T = fromSIunit_u(u3_rhoT(1 / (v3_ph(p4_T(T), h4V_p(p4_T(T)))), T));
                }
            }
            else
            {
                tempuV_T = 9999999999;  // error // need to check how to handle this
            }
            return tempuV_T;
        }

        public double uL_T(double T)
        {
            double tempuL_T = 0;
            T = toSIunit_T(T);
            if (T > 273.15 && T < 647.096)
            {
                if (T <= 623.15)
                {
                    tempuL_T = fromSIunit_u(u1_pT(p4_T(T), T));
                }
                else
                {
                    tempuL_T = fromSIunit_u(u3_rhoT(1 / (v3_ph(p4_T(T), h4L_p(p4_T(T)))), T));
                }
            }
            else
            {
                tempuL_T = 9999999999;  // error // need to check how to handle this
            }
            return tempuL_T;
        }

        public double u_pT(double p, double T)
        {
            double tempu_pT = 0;
            p = toSIunit_p(p);
            T = toSIunit_T(T);
            switch (region_pT(p, T))
            {
                case 1:
                    tempu_pT = fromSIunit_u(u1_pT(p, T));
                    break;
                case 2:
                    tempu_pT = fromSIunit_u(u2_pT(p, T));
                    break;
                case 3:
                    tempu_pT = fromSIunit_u(u3_rhoT(1 / v3_ph(p, h3_pT(p, T)), T));
                    break;
                case 4:
                    tempu_pT = 9999999999;  // error // need to check how to handle this
                    break;
                case 5:
                    tempu_pT = fromSIunit_u(u5_pT(p, T));
                    break;
                default:
                    tempu_pT = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempu_pT;
        }

        public double u_ph(double p, double h)
        {
            double tempu_ph = 0;
            double Ts = 0;
            double xs = 0;
            double u4v = 0;
            double u4L = 0;
            double v4V = 0;
            double v4L = 0;
            p = toSIunit_p(p);
            h = toSIunit_h(h);
            switch (region_ph(p, h))
            {
                case 1:
                    tempu_ph = fromSIunit_u(u1_pT(p, T1_ph(p, h)));
                    break;
                case 2:
                    tempu_ph = fromSIunit_u(u2_pT(p, T2_ph(p, h)));
                    break;
                case 3:
                    tempu_ph = fromSIunit_u(u3_rhoT(1 / v3_ph(p, h), T3_ph(p, h)));
                    break;
                case 4:
                    Ts = T4_p(p);
                    xs = x4_ph(p, h);
                    if (p < 16.529)
                    {
                        u4v = u2_pT(p, Ts);
                        u4L = u1_pT(p, Ts);
                    }
                    else
                    {
                        v4V = v3_ph(p, h4V_p(p));
                        u4v = u3_rhoT(1 / v4V, Ts);
                        v4L = v3_ph(p, h4L_p(p));
                        u4L = u3_rhoT(1 / v4L, Ts);
                    }
                    tempu_ph = fromSIunit_u((xs * u4v + (1 - xs) * u4L));
                    break;
                case 5:
                    Ts = T5_ph(p, h);
                    tempu_ph = fromSIunit_u(u5_pT(p, Ts));
                    break;
                default:
                    tempu_ph = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempu_ph;
        }

        public double u_ps(double p, double s)
        {
            double tempu_ps = 0;
            double x = 0;
        //    double u4v = 0;
            double uLp = 0;
            double uVp = 0;
         //   double u4L = 0;
            p = toSIunit_p(p);
            s = toSIunit_s(s);
            switch (region_ps(p, s))
            {
                case 1:
                    tempu_ps = fromSIunit_u(u1_pT(p, T1_ps(p, s)));
                    break;
                case 2:
                    tempu_ps = fromSIunit_u(u2_pT(p, T2_ps(p, s)));
                    break;
                case 3:
                    tempu_ps = fromSIunit_u(u3_rhoT(1 / v3_ps(p, s), T3_ps(p, s)));
                    break;
                case 4:
                    if (p < 16.529)
                    {
                        uLp = u1_pT(p, T4_p(p));
                        uVp = u2_pT(p, T4_p(p));
                    }
                    else
                    {
                        uLp = u3_rhoT(1 / (v3_ph(p, h4L_p(p))), T4_p(p));
                        uVp = u3_rhoT(1 / (v3_ph(p, h4V_p(p))), T4_p(p));
                    }
                    x = x4_ps(p, s);
                    tempu_ps = fromSIunit_u((x * uVp + (1 - x) * uLp));
                    break;
                case 5:
                    tempu_ps = fromSIunit_u(u5_pT(p, T5_ps(p, s)));
                    break;
                default:
                    tempu_ps = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempu_ps;
        }

        //***********************************************************************************************************
        //*1.9 Specific isobaric heat capacity (Cp)

        public double CpV_p(double p)
        {
            double tempCpV_p = 0;
            p = toSIunit_p(p);
            if (p > 0.000611657 && p < 22.06395)
            {
                if (p < 16.529)
                {
                    tempCpV_p = fromSIunit_Cp(Cp2_pT(p, T4_p(p)));
                }
                else
                {
                    tempCpV_p = fromSIunit_Cp(Cp3_rhoT(1 / (v3_ph(p, h4V_p(p))), T4_p(p)));
                }
            }
            else
            {
                tempCpV_p = 9999999999;  // error // need to check how to handle this
            }
            return tempCpV_p;
        }

        public double CpL_p(double p)
        {
            double tempCpL_p = 0;
            double T = 0;
            double h = 0;
            double v = 0;
            p = toSIunit_p(p);
            if (p > 0.000611657 && p < 22.06395)
            {
                if (p < 16.529)
                {
                    tempCpL_p = fromSIunit_Cp(Cp1_pT(p, T4_p(p)));
                }
                else
                {
                    T = T4_p(p);
                    h = h4L_p(p);
                    v = v3_ph(p, h4L_p(p));

                    tempCpL_p = fromSIunit_Cp(Cp3_rhoT(1 / (v3_ph(p, h4L_p(p))), T4_p(p)));
                }
            }
            else
            {
                tempCpL_p = 9999999999;  // error // need to check how to handle this
            }
            return tempCpL_p;
        }

        public double CpV_T(double T)
        {
            double tempCpV_T = 0;
            T = toSIunit_T(T);
            if (T > 273.15 && T < 647.096)
            {
                if (T <= 623.15)
                {
                    tempCpV_T = fromSIunit_Cp(Cp2_pT(p4_T(T), T));
                }
                else
                {
                    tempCpV_T = fromSIunit_Cp(Cp3_rhoT(1 / (v3_ph(p4_T(T), h4V_p(p4_T(T)))), T));
                }
            }
            else
            {
                tempCpV_T = 9999999999;  // error // need to check how to handle this
            }
            return tempCpV_T;
        }

        public double CpL_T(double T)
        {
            double tempCpL_T = 0;
            T = toSIunit_T(T);
            if (T > 273.15 && T < 647.096)
            {
                if (T <= 623.15)
                {
                    tempCpL_T = fromSIunit_Cp(Cp1_pT(p4_T(T), T));
                }
                else
                {
                    tempCpL_T = fromSIunit_Cp(Cp3_rhoT(1 / (v3_ph(p4_T(T), h4L_p(p4_T(T)))), T));
                }
            }
            else
            {
                tempCpL_T = 9999999999;  // error // need to check how to handle this
            }
            return tempCpL_T;
        }

        public double Cp_pT(double p, double T)
        {
            double tempCp_pT = 0;
            p = toSIunit_p(p);
            T = toSIunit_T(T);
            switch (region_pT(p, T))
            {
                case 1:
                    tempCp_pT = fromSIunit_Cp(Cp1_pT(p, T));
                    break;
                case 2:
                    tempCp_pT = fromSIunit_Cp(Cp2_pT(p, T));
                    break;
                case 3:
                    tempCp_pT = fromSIunit_Cp(Cp3_rhoT(1 / v3_ph(p, h3_pT(p, T)), T));
                    break;
                case 4:
                    tempCp_pT = 9999999999;  // error // need to check how to handle this
                    break;
                case 5:
                    tempCp_pT = fromSIunit_Cp(Cp5_pT(p, T));
                    break;
                default:
                    tempCp_pT = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempCp_pT;
        }

        public double Cp_ph(double p, double h)
        {
            double tempCp_ph = 0;
            p = toSIunit_p(p);
            h = toSIunit_h(h);
            switch (region_ph(p, h))
            {
                case 1:
                    tempCp_ph = fromSIunit_Cp(Cp1_pT(p, T1_ph(p, h)));
                    break;
                case 2:
                    tempCp_ph = fromSIunit_Cp(Cp2_pT(p, T2_ph(p, h)));
                    break;
                case 3:
                    tempCp_ph = fromSIunit_Cp(Cp3_rhoT(1 / v3_ph(p, h), T3_ph(p, h)));
                    break;
                case 4:
                    tempCp_ph = 9999999999;  // error // need to check how to handle this //#Not def. for mixture"
                    break;
                case 5:
                    tempCp_ph = fromSIunit_Cp(Cp5_pT(p, T5_ph(p, h)));
                    break;
                default:
                    tempCp_ph = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempCp_ph;
        }
        public double Cp_ps(double p, double s)
        {
            double tempCp_ps = 0;
            p = toSIunit_p(p);
            s = toSIunit_s(s);
            switch (region_ps(p, s))
            {
                case 1:
                    tempCp_ps = fromSIunit_Cp(Cp1_pT(p, T1_ps(p, s)));
                    break;
                case 2:
                    tempCp_ps = fromSIunit_Cp(Cp2_pT(p, T2_ps(p, s)));
                    break;
                case 3:
                    tempCp_ps = fromSIunit_Cp(Cp3_rhoT(1 / v3_ps(p, s), T3_ps(p, s)));
                    break;
                case 4:
                    tempCp_ps = 9999999999;  // error // need to check how to handle this //#Not def. for mixture"
                    break;
                case 5:
                    tempCp_ps = fromSIunit_Cp(Cp5_pT(p, T5_ps(p, s)));
                    break;
                default:
                    tempCp_ps = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempCp_ps;
        }
        //***********************************************************************************************************
        //*1.10 Specific isochoric heat capacity (Cv)

        public double CvV_p(double p)
        {
            double tempCvV_p = 0;
            p = toSIunit_p(p);
            if (p > 0.000611657 && p < 22.06395)
            {
                if (p < 16.529)
                {
                    tempCvV_p = fromSIunit_Cv(Cv2_pT(p, T4_p(p)));
                }
                else
                {
                    tempCvV_p = fromSIunit_Cv(Cv3_rhoT(1 / (v3_ph(p, h4V_p(p))), T4_p(p)));
                }
            }
            else
            {
                tempCvV_p = 9999999999;  // error // need to check how to handle this
            }
            return tempCvV_p;
        }

        public double CvL_p(double p)
        {
            double tempCvL_p = 0;
            p = toSIunit_p(p);
            if (p > 0.000611657 && p < 22.06395)
            {
                if (p < 16.529)
                {
                    tempCvL_p = fromSIunit_Cv(Cv1_pT(p, T4_p(p)));
                }
                else
                {
                    tempCvL_p = fromSIunit_Cv(Cv3_rhoT(1 / (v3_ph(p, h4L_p(p))), T4_p(p)));
                }
            }
            else
            {
                tempCvL_p = 9999999999;  // error // need to check how to handle this
            }
            return tempCvL_p;
        }

        public double CvV_T(double T)
        {
            double tempCvV_T = 0;
            T = toSIunit_T(T);
            if (T > 273.15 && T < 647.096)
            {
                if (T <= 623.15)
                {
                    tempCvV_T = fromSIunit_Cv(Cv2_pT(p4_T(T), T));
                }
                else
                {
                    tempCvV_T = fromSIunit_Cv(Cv3_rhoT(1 / (v3_ph(p4_T(T), h4V_p(p4_T(T)))), T));
                }
            }
            else
            {
                tempCvV_T = 9999999999;  // error // need to check how to handle this
            }
            return tempCvV_T;
        }

        public double CvL_T(double T)
        {
            double tempCvL_T = 0;
            T = toSIunit_T(T);
            if (T > 273.15 && T < 647.096)
            {
                if (T <= 623.15)
                {
                    tempCvL_T = fromSIunit_Cv(Cv1_pT(p4_T(T), T));
                }
                else
                {
                    tempCvL_T = fromSIunit_Cv(Cv3_rhoT(1 / (v3_ph(p4_T(T), h4L_p(p4_T(T)))), T));
                }
            }
            else
            {
                tempCvL_T = 9999999999;  // error // need to check how to handle this
            }
            return tempCvL_T;
        }

        public double Cv_pT(double p, double T)
        {
            double tempCv_pT = 0;
            p = toSIunit_p(p);
            T = toSIunit_T(T);
            switch (region_pT(p, T))
            {
                case 1:
                    tempCv_pT = fromSIunit_Cv(Cv1_pT(p, T));
                    break;
                case 2:
                    tempCv_pT = fromSIunit_Cv(Cv2_pT(p, T));
                    break;
                case 3:
                    tempCv_pT = fromSIunit_Cv(Cv3_rhoT(1 / v3_ph(p, h3_pT(p, T)), T));
                    break;
                case 4:
                    tempCv_pT = 9999999999;  // error // need to check how to handle this
                    break;
                case 5:
                    tempCv_pT = fromSIunit_Cv(Cv5_pT(p, T));
                    break;
                default:
                    tempCv_pT = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempCv_pT;
        }

        public double Cv_ph(double p, double h)
        {
            double tempCv_ph = 0;
            p = toSIunit_p(p);
            h = toSIunit_h(h);
            switch (region_ph(p, h))
            {
                case 1:
                    tempCv_ph = fromSIunit_Cv(Cv1_pT(p, T1_ph(p, h)));
                    break;
                case 2:
                    tempCv_ph = fromSIunit_Cv(Cv2_pT(p, T2_ph(p, h)));
                    break;
                case 3:
                    tempCv_ph = fromSIunit_Cv(Cv3_rhoT(1 / v3_ph(p, h), T3_ph(p, h)));
                    break;
                case 4:
                    tempCv_ph = 9999999999;  // error // need to check how to handle this //#Not def. for mixture"
                    break;
                case 5:
                    tempCv_ph = fromSIunit_Cv(Cv5_pT(p, T5_ph(p, h)));
                    break;
                default:
                    tempCv_ph = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempCv_ph;
        }

        public double Cv_ps(double p, double s)
        {
            double tempCv_ps = 0;
            p = toSIunit_p(p);
            s = toSIunit_s(s);
            switch (region_ps(p, s))
            {
                case 1:
                    tempCv_ps = fromSIunit_Cv(Cv1_pT(p, T1_ps(p, s)));
                    break;
                case 2:
                    tempCv_ps = fromSIunit_Cv(Cv2_pT(p, T2_ps(p, s)));
                    break;
                case 3:
                    tempCv_ps = fromSIunit_Cv(Cv3_rhoT(1 / v3_ps(p, s), T3_ps(p, s)));
                    break;
                case 4:
                    tempCv_ps = 9999999999;  // error // need to check how to handle this //#Not def. for mixture
                    break;
                case 5:
                    tempCv_ps = fromSIunit_Cv(Cv5_pT(p, T5_ps(p, s)));
                    break;
                default:
                    tempCv_ps = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempCv_ps;
        }

        //***********************************************************************************************************
        //*1.11 Speed of sound


        public double wV_p(double p)
        {
            double tempwV_p = 0;
            p = toSIunit_p(p);
            if (p > 0.000611657 && p < 22.06395)
            {
                if (p < 16.529)
                {
                    tempwV_p = fromSIunit_w(w2_pT(p, T4_p(p)));
                }
                else
                {
                    tempwV_p = fromSIunit_w(w3_rhoT(1 / (v3_ph(p, h4V_p(p))), T4_p(p)));
                }
            }
            else
            {
                tempwV_p = 9999999999;  // error // need to check how to handle this
            }
            return tempwV_p;
        }

        public double wL_p(double p)
        {
            double tempwL_p = 0;
            p = toSIunit_p(p);
            if (p > 0.000611657 && p < 22.06395)
            {
                if (p < 16.529)
                {
                    tempwL_p = fromSIunit_w(w1_pT(p, T4_p(p)));
                }
                else
                {
                    tempwL_p = fromSIunit_w(w3_rhoT(1 / (v3_ph(p, h4L_p(p))), T4_p(p)));
                }
            }
            else
            {
                tempwL_p = 9999999999;  // error // need to check how to handle this
            }
            return tempwL_p;
        }

        public double wV_T(double T)
        {
            double tempwV_T = 0;
            T = toSIunit_T(T);
            if (T > 273.15 && T < 647.096)
            {
                if (T <= 623.15)
                {
                    tempwV_T = fromSIunit_w(w2_pT(p4_T(T), T));
                }
                else
                {
                    tempwV_T = fromSIunit_w(w3_rhoT(1 / (v3_ph(p4_T(T), h4V_p(p4_T(T)))), T));
                }
            }
            else
            {
                tempwV_T = 9999999999;  // error // need to check how to handle this
            }
            return tempwV_T;
        }

        public double wL_T(double T)
        {
            double tempwL_T = 0;
            T = toSIunit_T(T);
            if (T > 273.15 && T < 647.096)
            {
                if (T <= 623.15)
                {
                    tempwL_T = fromSIunit_w(w1_pT(p4_T(T), T));
                }
                else
                {
                    tempwL_T = fromSIunit_w(w3_rhoT(1 / (v3_ph(p4_T(T), h4L_p(p4_T(T)))), T));
                }
            }
            else
            {
                tempwL_T = 9999999999;  // error // need to check how to handle this
            }
            return tempwL_T;
        }

        public double w_pT(double p, double T)
        {
            double tempw_pT = 0;
            p = toSIunit_p(p);
            T = toSIunit_T(T);
            switch (region_pT(p, T))
            {
                case 1:
                    tempw_pT = fromSIunit_w(w1_pT(p, T));
                    break;
                case 2:
                    tempw_pT = fromSIunit_w(w2_pT(p, T));
                    break;
                case 3:
                    tempw_pT = fromSIunit_w(w3_rhoT(1 / v3_ph(p, h3_pT(p, T)), T));
                    break;
                case 4:
                    tempw_pT = 9999999999;  // error // need to check how to handle this
                    break;
                case 5:
                    tempw_pT = fromSIunit_w(w5_pT(p, T));
                    break;
                default:
                    tempw_pT = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempw_pT;
        }

        public double w_ph(double p, double h)
        {
            double tempw_ph = 0;
            p = toSIunit_p(p);
            h = toSIunit_h(h);
            switch (region_ph(p, h))
            {
                case 1:
                    tempw_ph = fromSIunit_w(w1_pT(p, T1_ph(p, h)));
                    break;
                case 2:
                    tempw_ph = fromSIunit_w(w2_pT(p, T2_ph(p, h)));
                    break;
                case 3:
                    tempw_ph = fromSIunit_w(w3_rhoT(1 / v3_ph(p, h), T3_ph(p, h)));
                    break;
                case 4:
                    tempw_ph = 9999999999;  // error // need to check how to handle this //#Not def. for mixture
                    break;
                case 5:
                    tempw_ph = fromSIunit_w(w5_pT(p, T5_ph(p, h)));
                    break;
                default:
                    tempw_ph = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempw_ph;
        }

        public double w_ps(double p, double s)
        {
            double tempw_ps = 0;
            p = toSIunit_p(p);
            s = toSIunit_s(s);
            switch (region_ps(p, s))
            {
                case 1:
                    tempw_ps = fromSIunit_w(w1_pT(p, T1_ps(p, s)));
                    break;
                case 2:
                    tempw_ps = fromSIunit_w(w2_pT(p, T2_ps(p, s)));
                    break;
                case 3:
                    tempw_ps = fromSIunit_w(w3_rhoT(1 / v3_ps(p, s), T3_ps(p, s)));
                    break;
                case 4:
                    tempw_ps = 9999999999;  // error // need to check how to handle this //#Not def. for mixture
                    break;
                case 5:
                    tempw_ps = fromSIunit_w(w5_pT(p, T5_ps(p, s)));
                    break;
                default:
                    tempw_ps = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempw_ps;
        }

        //***********************************************************************************************************
        //*1.12 Viscosity

        public double my_pT(double p, double T)
        {
            double tempmy_pT = 0;
            p = toSIunit_p(p);
            T = toSIunit_T(T);
            switch (region_pT(p, T))
            {
                case 4:
                    tempmy_pT = 9999999999;  // error // need to check how to handle this
                    break;
                case 1:
                case 2:
                case 3:
                case 5:
                    tempmy_pT = fromSIunit_my(my_AllRegions_pT(p, T));
                    break;
                default:
                    tempmy_pT = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempmy_pT;
        }

        public double my_ph(double p, double h)
        {
            double tempmy_ph = 0;
            p = toSIunit_p(p);
            h = toSIunit_h(h);
            switch (region_ph(p, h))
            {
                case 1:
                case 2:
                case 3:
                case 5:
                    tempmy_ph = fromSIunit_my(my_AllRegions_ph(p, h));
                    break;
                case 4:
                    tempmy_ph = 9999999999;  // error // need to check how to handle this
                    break;
                default:
                    tempmy_ph = 9999999999;  // error // need to check how to handle this
                    break;
            }
            return tempmy_ph;
        }

        public double my_ps(double p, double s)
        {
            return my_ph(p, h_ps(p, s));
        }

        public double Pr_pT(double p, double T)
        {
            double Cp = 0;
            double my = 0;
            double tc = 0;
            Cp = toSIunit_Cp(Cp_pT(p, T));
            my = toSIunit_my(my_pT(p, T));
            tc = toSIunit_tc(tc_pT(p, T));
            return Cp * 1000 * my / tc;
        }

        public double Pr_ph(double p, double h)
        {
            double Cp = 0;
            double my = 0;
            double tc = 0;
            Cp = toSIunit_Cp(Cp_ph(p, h));
            my = toSIunit_my(my_ph(p, h));
            tc = toSIunit_tc(tc_ph(p, h));
            return Cp * 1000 * my / tc;
        }

        public double Kappa_pT(double p, double T)
        {
            double Cp = 0;
            double Cv = 0;
            Cp = Cp_pT(p, T);
            Cv = Cv_pT(p, T);
            return Cp / Cv;
        }

        public double Kappa_ph(double p, double h)
        {
            double Cp = 0;
            double Cv;
            Cv = Cv_ph(p, h);
            Cp = Cp_ph(p, h);
            return Cp / Cv;
        }
        //***********************************************************************************************************
        //*1.15 Surface tension

        public double st_t(double T)
        {
            T = toSIunit_T(T);
            return fromSIunit_st(Surface_Tension_T(T));
        }

        public double st_p(double p)
        {
            double T;
            T = Tsat_p(p);
            T = toSIunit_T(T);
            return fromSIunit_st(Surface_Tension_T(T));
        }

        //***********************************************************************************************************
        //*1.16 Thermal conductivity
        public double tcL_p(double p)
        {
            double T = 0;
            double v = 0;
            T = Tsat_p(p);
            v = vL_p(p);
            p = toSIunit_p(p);
            T = toSIunit_T(T);
            v = toSIunit_v(v);
            return fromSIunit_tc(tc_ptrho(p, T, 1 / v));
        }

        public double tcV_p(double p)
        {
            double T = 0;
            double v = 0;
            T = Tsat_p(p);
            v = vV_p(p);
            p = toSIunit_p(p);
            T = toSIunit_T(T);
            v = toSIunit_v(v);
            return fromSIunit_tc(tc_ptrho(p, T, 1 / v));
        }

        public double tcL_T(double T)
        {
            double p = 0;
            double v = 0;
            p = psat_T(T);
            v = vL_T(T);
            p = toSIunit_p(p);
            T = toSIunit_T(T);
            v = toSIunit_v(v);
            return fromSIunit_tc(tc_ptrho(p, T, 1 / v));
        }

        public double tcV_T(double T)
        {
            double p = 0;
            double v = 0;
            p = psat_T(T);
            v = vV_T(T);
            p = toSIunit_p(p);
            T = toSIunit_T(T);
            v = toSIunit_v(v);
            return fromSIunit_tc(tc_ptrho(p, T, 1 / v));
        }

        public double tc_pT(double p, double T)
        {
            double v;
            v = v_pT(p, T);
            p = toSIunit_p(p);
            T = toSIunit_T(T);
            v = toSIunit_v(v);
            return fromSIunit_tc(tc_ptrho(p, T, 1 / v));
        }

        public double tc_ph(double p, double h)
        {
            double v = 0;
            double T = 0;
            v = v_ph(p, h);
            T = T_ph(p, h);
            p = toSIunit_p(p);
            T = toSIunit_T(T);
            v = toSIunit_v(v);
            return fromSIunit_tc(tc_ptrho(p, T, 1 / v));
        }

        public double tc_hs(double h, double s)
        {
            double p = 0;
            double v = 0;
            double T = 0;
            p = p_hs(h, s);
            v = v_ph(p, h);
            T = T_ph(p, h);
            p = toSIunit_p(p);
            T = toSIunit_T(T);
            v = toSIunit_v(v);
            return fromSIunit_tc(tc_ptrho(p, T, 1 / v));
        }

        //***********************************************************************************************************
        //*1.17 Vapour fraction

        public double x_ph(double p, double h)
        {
            double tempx_ph = 0;
            p = toSIunit_p(p);
            h = toSIunit_h(h);
            if (p > 0.000611657 && p < 22.06395)
            {
                tempx_ph = fromSIunit_x(x4_ph(p, h));
            }
            else
            {
                tempx_ph = 9999999999;  // error // need to check how to handle this
            }
            return tempx_ph;
        }

        public double x_ps(double p, double s)
        {
            double tempx_ps = 0;
            p = toSIunit_p(p);
            s = toSIunit_s(s);
            if (p > 0.000611657 && p < 22.06395)
            {
                tempx_ps = fromSIunit_x(x4_ps(p, s));
            }
            else
            {
                tempx_ps = 9999999999;  // error // need to check how to handle this
            }
            return tempx_ps;
        }

        //***********************************************************************************************************
        //*1.18 Vapour Volume Fraction

        public double vx_ph(double p, double h)
        {
            double tempvx_ph = 0;
            double vL = 0;
            double vV = 0;
            double xs = 0;
            p = toSIunit_p(p);
            h = toSIunit_h(h);
            if (p > 0.000611657 && p < 22.06395)
            {
                if (p < 16.529)
                {
                    vL = v1_pT(p, T4_p(p));
                    vV = v2_pT(p, T4_p(p));
                }
                else
                {
                    vL = v3_ph(p, h4L_p(p));
                    vV = v3_ph(p, h4V_p(p));
                }
                xs = x4_ph(p, h);
                tempvx_ph = fromSIunit_vx((xs * vV / (xs * vV + (1 - xs) * vL)));
            }
            else
            {
                tempvx_ph = 9999999999;  // error // need to check how to handle this
            }
            return tempvx_ph;
        }

        public double vx_ps(double p, double s)
        {
            double tempvx_ps = 0;
            double vL = 0;
            double vV = 0;
            double xs = 0;
            p = toSIunit_p(p);
            s = toSIunit_s(s);
            if (p > 0.000611657 && p < 22.06395)
            {
                if (p < 16.529)
                {
                    vL = v1_pT(p, T4_p(p));
                    vV = v2_pT(p, T4_p(p));
                }
                else
                {
                    vL = v3_ph(p, h4L_p(p));
                    vV = v3_ph(p, h4V_p(p));
                }
                xs = x4_ps(p, s);
                tempvx_ps = fromSIunit_vx((xs * vV / (xs * vV + (1 - xs) * vL)));
            }
            else
            {
                tempvx_ps = 9999999999;  // error // need to check how to handle this
            }
            return tempvx_ps;
        }

        //***********************************************************************************************************
        //*2 IAPWS IF 97 Calling functions                                                                          *
        //***********************************************************************************************************
        //
        //***********************************************************************************************************
        //*2.1 Functions for region 1

        public double v1_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //5 Equations for Region 1, Section. 5.1 Basic Equation
            //Eqution 7, Table 3, Page 6
            int i = 0;
            double ps = 0;
            double tau = 0;
            double g_p = 0;
            //Variant *I1 = nullptr;
            int[] I1 = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32 };
            int[] J1 = { -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41 };
            double[] n1 = { 0.14632971213167, -0.84548187169114, -3.756360367204, 3.3855169168385, -0.95791963387872, 0.15772038513228, -0.016616417199501, 8.1214629983568E-04, 2.8319080123804E-04, -6.0706301565874E-04, -0.018990068218419, -0.032529748770505, -0.021841717175414, -5.283835796993E-05, -4.7184321073267E-04, -3.0001780793026E-04, 4.7661393906987E-05, -4.4141845330846E-06, -7.2694996297594E-16, -3.1679644845054E-05, -2.8270797985312E-06, -8.5205128120103E-10, -2.2425281908E-06, -6.5171222895601E-07, -1.4341729937924E-13, -4.0516996860117E-07, -1.2734301741641E-09, -1.7424871230634E-10, -6.8762131295531E-19, 1.4478307828521E-20, 2.6335781662795E-23, -1.1947622640071E-23, 1.8228094581404E-24, -9.3537087292458E-26 };
            double R = 0.461526; //kJ/(kg K)

            // I1[34] = 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32;
            //J1[34] = -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41;
            //  n1[34]= 0.14632971213167, -0.84548187169114, -3.756360367204, 3.3855169168385, -0.95791963387872, 0.15772038513228, -0.016616417199501, 8.1214629983568E-04, 2.8319080123804E-04, -6.0706301565874E-04, -0.018990068218419, -0.032529748770505, -0.021841717175414, -5.283835796993E-05, -4.7184321073267E-04, -3.0001780793026E-04, 4.7661393906987E-05, -4.4141845330846E-06, -7.2694996297594E-16, -3.1679644845054E-05, -2.8270797985312E-06, -8.5205128120103E-10, -2.2425281908E-06, -6.5171222895601E-07, -1.4341729937924E-13, -4.0516996860117E-07, -1.2734301741641E-09, -1.7424871230634E-10, -6.8762131295531E-19, 1.4478307828521E-20, 2.6335781662795E-23, -1.1947622640071E-23, 1.8228094581404E-24, -9.3537087292458E-26;
            ps = p / 16.53;
            tau = 1386 / T;
            g_p = 0;
            for (i = 0; i <= 33; i++)
            {
                g_p = g_p - n1[i] * I1[i] * Math.Pow(7.1 - ps, I1[i] - 1) * Math.Pow(tau - 1.222, J1[i]);
            }
            return R * T / p * ps * g_p / 1000;
        }

        public double h1_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //5 Equations for Region 1, Section. 5.1 Basic Equation
            //Eqution 7, Table 3, Page 6
            int i = 0;
        //    double ps = 0;
            double tau = 0;
            double g_t = 0;
            int[] I1 = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32 };
            int[] J1 = { -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41 };
            double[] n1 = { 0.14632971213167, -0.84548187169114, -3.756360367204, 3.3855169168385, -0.95791963387872, 0.15772038513228, -0.016616417199501, 8.1214629983568E-04, 2.8319080123804E-04, -6.0706301565874E-04, -0.018990068218419, -0.032529748770505, -0.021841717175414, -5.283835796993E-05, -4.7184321073267E-04, -3.0001780793026E-04, 4.7661393906987E-05, -4.4141845330846E-06, -7.2694996297594E-16, -3.1679644845054E-05, -2.8270797985312E-06, -8.5205128120103E-10, -2.2425281908E-06, -6.5171222895601E-07, -1.4341729937924E-13, -4.0516996860117E-07, -1.2734301741641E-09, -1.7424871230634E-10, -6.8762131295531E-19, 1.4478307828521E-20, 2.6335781662795E-23, -1.1947622640071E-23, 1.8228094581404E-24, -9.3537087292458E-26 };

            double R = 0.461526; //kJ/(kg K)
            p = p / 16.53;
            tau = 1386 / T;
            g_t = 0;
            for (i = 0; i <= 33; i++)
            {
                g_t = g_t + (n1[i] * Math.Pow(7.1 - p, I1[i]) * J1[i] * Math.Pow(tau - 1.222, J1[i] - 1));
            }
            return R * T * tau * g_t;
        }

        public double u1_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //5 Equations for Region 1, Section. 5.1 Basic Equation
            //Eqution 7, Table 3, Page 6
            int i = 0;
            double tau = 0;
            double g_t = 0;
            double g_p = 0;
            int[] I1 = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32 };
            int[] J1 = { -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41 };
            double[] n1 = { 0.14632971213167, -0.84548187169114, -3.756360367204, 3.3855169168385, -0.95791963387872, 0.15772038513228, -0.016616417199501, 8.1214629983568E-04, 2.8319080123804E-04, -6.0706301565874E-04, -0.018990068218419, -0.032529748770505, -0.021841717175414, -5.283835796993E-05, -4.7184321073267E-04, -3.0001780793026E-04, 4.7661393906987E-05, -4.4141845330846E-06, -7.2694996297594E-16, -3.1679644845054E-05, -2.8270797985312E-06, -8.5205128120103E-10, -2.2425281908E-06, -6.5171222895601E-07, -1.4341729937924E-13, -4.0516996860117E-07, -1.2734301741641E-09, -1.7424871230634E-10, -6.8762131295531E-19, 1.4478307828521E-20, 2.6335781662795E-23, -1.1947622640071E-23, 1.8228094581404E-24, -9.3537087292458E-26 };
            double R = 0.461526; //kJ/(kg K)

            p = p / 16.53;
            tau = 1386 / T;
            g_t = 0;
            g_p = 0;
            for (i = 0; i <= 33; i++)
            {
                g_p = g_p - n1[i] * I1[i] * Math.Pow(7.1 - p, I1[i] - 1) * Math.Pow(tau - 1.222, J1[i]);
                g_t = g_t + (n1[i] * Math.Pow(7.1 - p, I1[i]) * J1[i] * Math.Pow(tau - 1.222, J1[i] - 1));
            }
            return R * T * (tau * g_t - p * g_p);
        }

        public double s1_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //5 Equations for Region 1, Section. 5.1 Basic Equation
            //Eqution 7, Table 3, Page 6
            int i = 0;
            double g = 0;
            double g_t = 0;
            double R = 0.461526; //kJ/(kg K)
            int[] I1 = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32 };
            int[] J1 = { -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41 };
            double[] n1 = { 0.14632971213167, -0.84548187169114, -3.756360367204, 3.3855169168385, -0.95791963387872, 0.15772038513228, -0.016616417199501, 8.1214629983568E-04, 2.8319080123804E-04, -6.0706301565874E-04, -0.018990068218419, -0.032529748770505, -0.021841717175414, -5.283835796993E-05, -4.7184321073267E-04, -3.0001780793026E-04, 4.7661393906987E-05, -4.4141845330846E-06, -7.2694996297594E-16, -3.1679644845054E-05, -2.8270797985312E-06, -8.5205128120103E-10, -2.2425281908E-06, -6.5171222895601E-07, -1.4341729937924E-13, -4.0516996860117E-07, -1.2734301741641E-09, -1.7424871230634E-10, -6.8762131295531E-19, 1.4478307828521E-20, 2.6335781662795E-23, -1.1947622640071E-23, 1.8228094581404E-24, -9.3537087292458E-26 };

            p = p / 16.53;
            T = 1386 / T;
            g = 0;
            g_t = 0;
            for (i = 0; i <= 33; i++)
            {
                g_t = g_t + (n1[i] * Math.Pow(7.1 - p, I1[i]) * J1[i] * Math.Pow(T - 1.222, J1[i] - 1));
                g = g + n1[i] * Math.Pow(7.1 - p, I1[i]) * Math.Pow(T - 1.222, J1[i]);
            }
            return R * T * g_t - R * g;
        }

        public double Cp1_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //5 Equations for Region 1, Section. 5.1 Basic Equation
            //Eqution 7, Table 3, Page 6
            int i = 0;
            double G_tt = 0;
            int[] I1 = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32 };
            int[] J1 = { -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41 };
            double[] n1 = { 0.14632971213167, -0.84548187169114, -3.756360367204, 3.3855169168385, -0.95791963387872, 0.15772038513228, -0.016616417199501, 8.1214629983568E-04, 2.8319080123804E-04, -6.0706301565874E-04, -0.018990068218419, -0.032529748770505, -0.021841717175414, -5.283835796993E-05, -4.7184321073267E-04, -3.0001780793026E-04, 4.7661393906987E-05, -4.4141845330846E-06, -7.2694996297594E-16, -3.1679644845054E-05, -2.8270797985312E-06, -8.5205128120103E-10, -2.2425281908E-06, -6.5171222895601E-07, -1.4341729937924E-13, -4.0516996860117E-07, -1.2734301741641E-09, -1.7424871230634E-10, -6.8762131295531E-19, 1.4478307828521E-20, 2.6335781662795E-23, -1.1947622640071E-23, 1.8228094581404E-24, -9.3537087292458E-26 };
            double R = 0.461526; //kJ/(kg K)

            p = p / 16.53;
            T = 1386 / T;
            G_tt = 0;
            for (i = 0; i <= 33; i++)
            {
                G_tt = G_tt + (n1[i] * Math.Pow(7.1 - p, I1[i]) * J1[i] * (J1[i] - 1) * Math.Pow(T - 1.222, J1[i] - 2));
            }
            return -R * Math.Pow(T, 2) * G_tt;

        }

        public double Cv1_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //5 Equations for Region 1, Section. 5.1 Basic Equation
            //Eqution 7, Table 3, Page 6
            int i = 0;
            double g_p = 0;
            double g_pp = 0;
            double g_pt = 0;
            double G_tt = 0;
            double R = 0.461526; //kJ/(kg K)
            int[] I1 = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32 };
            int[] J1 = { -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41 };
            double[] n1 = { 0.14632971213167, -0.84548187169114, -3.756360367204, 3.3855169168385, -0.95791963387872, 0.15772038513228, -0.016616417199501, 8.1214629983568E-04, 2.8319080123804E-04, -6.0706301565874E-04, -0.018990068218419, -0.032529748770505, -0.021841717175414, -5.283835796993E-05, -4.7184321073267E-04, -3.0001780793026E-04, 4.7661393906987E-05, -4.4141845330846E-06, -7.2694996297594E-16, -3.1679644845054E-05, -2.8270797985312E-06, -8.5205128120103E-10, -2.2425281908E-06, -6.5171222895601E-07, -1.4341729937924E-13, -4.0516996860117E-07, -1.2734301741641E-09, -1.7424871230634E-10, -6.8762131295531E-19, 1.4478307828521E-20, 2.6335781662795E-23, -1.1947622640071E-23, 1.8228094581404E-24, -9.3537087292458E-26 };
            p = p / 16.53;
            T = 1386 / T;
            g_p = 0;
            g_pp = 0;
            g_pt = 0;
            G_tt = 0;
            for (i = 0; i <= 33; i++)
            {
                g_p = g_p - n1[i] * I1[i] * Math.Pow(7.1 - p, I1[i] - 1) * Math.Pow(T - 1.222, J1[i]);
                g_pp = g_pp + n1[i] * I1[i] * (I1[i] - 1) * Math.Pow(7.1 - p, I1[i] - 2) * Math.Pow(T - 1.222, J1[i]);
                g_pt = g_pt - n1[i] * I1[i] * Math.Pow(7.1 - p, I1[i] - 1) * J1[i] * Math.Pow(T - 1.222, J1[i] - 1);
                G_tt = G_tt + n1[i] * Math.Pow(7.1 - p, I1[i]) * J1[i] * (J1[i] - 1) * Math.Pow(T - 1.222, J1[i] - 2);
            }
            return R * (-(Math.Pow(T, 2) * G_tt) + Math.Pow(g_p - T * g_pt, 2) / g_pp);
        }

        public double w1_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //5 Equations for Region 1, Section. 5.1 Basic Equation
            //Eqution 7, Table 3, Page 6
            int i = 0;
            double g_p = 0;
            double g_pp = 0;
            double g_pt = 0;
            double G_tt = 0;
            double tau = 0;
            int[] I1 = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32 };
            int[] J1 = { -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41 };
            double[] n1 = { 0.14632971213167, -0.84548187169114, -3.756360367204, 3.3855169168385, -0.95791963387872, 0.15772038513228, -0.016616417199501, 8.1214629983568E-04, 2.8319080123804E-04, -6.0706301565874E-04, -0.018990068218419, -0.032529748770505, -0.021841717175414, -5.283835796993E-05, -4.7184321073267E-04, -3.0001780793026E-04, 4.7661393906987E-05, -4.4141845330846E-06, -7.2694996297594E-16, -3.1679644845054E-05, -2.8270797985312E-06, -8.5205128120103E-10, -2.2425281908E-06, -6.5171222895601E-07, -1.4341729937924E-13, -4.0516996860117E-07, -1.2734301741641E-09, -1.7424871230634E-10, -6.8762131295531E-19, 1.4478307828521E-20, 2.6335781662795E-23, -1.1947622640071E-23, 1.8228094581404E-24, -9.3537087292458E-26 };
            double R = 0.461526; //kJ/(kg K)
            p = p / 16.53;
            tau = 1386 / T;
            g_p = 0;
            g_pp = 0;
            g_pt = 0;
            G_tt = 0;
            for (i = 0; i <= 33; i++)
            {
                g_p = g_p - n1[i] * I1[i] * Math.Pow(7.1 - p, I1[i] - 1) * Math.Pow(tau - 1.222, J1[i]);
                g_pp = g_pp + n1[i] * I1[i] * (I1[i] - 1) * Math.Pow(7.1 - p, I1[i] - 2) * Math.Pow(tau - 1.222, J1[i]);
                g_pt = g_pt - n1[i] * I1[i] * Math.Pow(7.1 - p, I1[i] - 1) * J1[i] * Math.Pow(tau - 1.222, J1[i] - 1);
                G_tt = G_tt + n1[i] * Math.Pow(7.1 - p, I1[i]) * J1[i] * (J1[i] - 1) * Math.Pow(tau - 1.222, J1[i] - 2);
            }
            return Math.Pow(1000 * R * T * Math.Pow(g_p, 2) / (Math.Pow(g_p - tau * g_pt, 2) / (Math.Pow(tau, 2) * G_tt) - g_pp), 0.5);
        }

        public double T1_ph(double p, double h)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //5 Equations for Region 1, Section. 5.1 Basic Equation, 5.2.1 The Backward Equation T ( p,h )
            //Eqution 11, Table 6, Page 10
            int i = 0;
            double T = 0;
            int[] I1 = { 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6 };
            int[] J1 = { 0, 1, 2, 6, 22, 32, 0, 1, 2, 3, 4, 10, 32, 10, 32, 10, 32, 32, 32, 32 };
            double[] n1 = { -238.72489924521, 404.21188637945, 113.49746881718, -5.8457616048039, -1.528548241314E-04, -1.0866707695377E-06, -13.391744872602, 43.211039183559, -54.010067170506, 30.535892203916, -6.5964749423638, 9.3965400878363E-03, 1.157364750534E-07, -2.5858641282073E-05, -4.0644363084799E-09, 6.6456186191635E-08, 8.0670734103027E-11, -9.3477771213947E-13, 5.8265442020601E-15, -1.5020185953503E-17 };

            h = h / 2500;
            T = 0;
            for (i = 0; i <= 19; i++)
            {
                T = T + n1[i] * Math.Pow(p, I1[i]) * Math.Pow(h + 1, J1[i]);
            }
            return T;
        }
        public double T1_ps(double p, double s)
        {
            double tempT1_ps = 0;
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //5 Equations for Region 1, Section. 5.1 Basic Equation, 5.2.2 The Backward Equation T ( p, s )
            //Eqution 13, Table 8, Page 11
            int i = 0;
            int[] I1 = { 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4 };
            int[] J1 = { 0, 1, 2, 3, 11, 31, 0, 1, 2, 3, 12, 31, 0, 1, 2, 9, 31, 10, 32, 32 };
            double[] n1 = { 174.78268058307, 34.806930892873, 6.5292584978455, 0.33039981775489, -1.9281382923196E-07, -2.4909197244573E-23, -0.26107636489332, 0.22592965981586, -0.064256463395226, 7.8876289270526E-03, 3.5672110607366E-10, 1.7332496994895E-24, 5.6608900654837E-04, -3.2635483139717E-04, 4.4778286690632E-05, -5.1322156908507E-10, -4.2522657042207E-26, 2.6400441360689E-13, 7.8124600459723E-29, -3.0732199903668E-31 };

            tempT1_ps = 0;
            for (i = 0; i <= 19; i++)
            {
                tempT1_ps = tempT1_ps + n1[i] * Math.Pow(p, I1[i]) * Math.Pow(s + 2, J1[i]);
            }
            return tempT1_ps;
        }

        public double p1_hs(double h, double s)
        {
            //Supplementary Release on Backward Equations for Pressure as a Function of Enthalpy and Entropy p(h,s) to the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam
            //5 Backward Equation p(h,s) for Region 1
            //Eqution 1, Table 2, Page 5
            int i = 0;
            double p = 0;
            int[] I1 = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 5 };
            int[] J1 = { 0, 1, 2, 4, 5, 6, 8, 14, 0, 1, 4, 6, 0, 1, 10, 4, 1, 4, 0 };
            double[] n1 = { -0.691997014660582, -18.361254878756, -9.28332409297335, 65.9639569909906, -16.2060388912024, 450.620017338667, 854.68067822417, 6075.23214001162, 32.6487682621856, -26.9408844582931, -319.9478483343, -928.35430704332, 30.3634537455249, -65.0540422444146, -4309.9131651613, -747.512324096068, 730.000345529245, 1142.84032569021, -436.407041874559 };

            h = h / 3400;
            s = s / 7.6;
            p = 0;
            for (i = 0; i <= 18; i++)
            {
                p = p + n1[i] * Math.Pow(h + 0.05, I1[i]) * Math.Pow(s + 0.05, J1[i]);
            }
            return p * 100;
        }

        public double T1_prho(double p, double rho)
        {
            //Solve by iteration. Observe that fo low temperatures this equation has 2 solutions.
            //Solve with half interval method
         //   int i = 0;
            double Ts = 0;
            double Low_Bound = 0;
            double High_Bound = 0;
            double rhos = 0;
            Low_Bound = 273.15;
            High_Bound = T4_p(p);
            while (Math.Abs(rho - rhos) > 0.00001)
            {
                Ts = (Low_Bound + High_Bound) / 2;
                rhos = 1 / v1_pT(p, Ts);
                if (rhos < rho)
                {
                    High_Bound = Ts;
                }
                else
                {
                    Low_Bound = Ts;
                }
            }
            return Ts;
        }



        //***********************************************************************************************************
        //*2.2 Functions for region 2

        public double v2_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //6 Equations for Region 2, Section. 6.1 Basic Equation
            //Table 11 and 12, Page 14 and 15
            int i = 0;
            double tau = 0;
            double g0_pi = 0;
            double gr_pi = 0;
            int[] Ir = { -1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24 };
            int[] Jr = { -0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58 };
            double[] nr = { -1.7731742473213E-03, -0.017834862292358, -0.045996013696365, -0.057581259083432, -0.05032527872793, -3.3032641670203E-05, -1.8948987516315E-04, -3.9392777243355E-03, -0.043797295650573, -2.6674547914087E-05, 2.0481737692309E-08, 4.3870667284435E-07, -3.227767723857E-05, -1.5033924542148E-03, -0.040668253562649, -7.8847309559367E-10, 1.2790717852285E-08, 4.8225372718507E-07, 2.2922076337661E-06, -1.6714766451061E-11, -2.1171472321355E-03, -23.895741934104, -5.905956432427E-18, -1.2621808899101E-06, -0.038946842435739, 1.1256211360459E-11, -8.2311340897998, 1.9809712802088E-08, 1.0406965210174E-19, -1.0234747095929E-13, -1.0018179379511E-09, -8.0882908646985E-11, 0.10693031879409, -0.33662250574171, 8.9185845355421E-25, 3.0629316876232E-13, -4.2002467698208E-06, -5.9056029685639E-26, 3.7826947613457E-06, -1.2768608934681E-15, 7.3087610595061E-29, 5.5414715350778E-17, -9.436970724121E-07 };
            int[] J0 = { 0, 1, -5, -4, -3, -2, -1, 2, 3 };
            double[] n0 = { -9.6927686500217, 10.086655968018, -0.005608791128302, 0.071452738081455, -0.40710498223928, 1.4240819171444, -4.383951131945, -0.28408632460772, 0.021268463753307 };
            double R = 0.461526; //kJ/(kg K)

            tau = 540 / T;
            g0_pi = 1 / p;
            gr_pi = 0;
            for (i = 0; i <= 42; i++)
            {
                gr_pi = gr_pi + nr[i] * Ir[i] * Math.Pow(p, Ir[i] - 1) * Math.Pow(tau - 0.5, Jr[i]);
            }
            return R * T / p * p * (g0_pi + gr_pi) / 1000;
        }

        public double h2_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //6 Equations for Region 2, Section. 6.1 Basic Equation
            //Table 11 and 12, Page 14 and 15
            int i = 0;
            double tau = 0;
            double g0_tau = 0;
            double gr_tau = 0;
            int[] Ir = { -1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24 };
            int[] Jr = { -0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58 };
            double[] nr = { -1.7731742473213E-03, -0.017834862292358, -0.045996013696365, -0.057581259083432, -0.05032527872793, -3.3032641670203E-05, -1.8948987516315E-04, -3.9392777243355E-03, -0.043797295650573, -2.6674547914087E-05, 2.0481737692309E-08, 4.3870667284435E-07, -3.227767723857E-05, -1.5033924542148E-03, -0.040668253562649, -7.8847309559367E-10, 1.2790717852285E-08, 4.8225372718507E-07, 2.2922076337661E-06, -1.6714766451061E-11, -2.1171472321355E-03, -23.895741934104, -5.905956432427E-18, -1.2621808899101E-06, -0.038946842435739, 1.1256211360459E-11, -8.2311340897998, 1.9809712802088E-08, 1.0406965210174E-19, -1.0234747095929E-13, -1.0018179379511E-09, -8.0882908646985E-11, 0.10693031879409, -0.33662250574171, 8.9185845355421E-25, 3.0629316876232E-13, -4.2002467698208E-06, -5.9056029685639E-26, 3.7826947613457E-06, -1.2768608934681E-15, 7.3087610595061E-29, 5.5414715350778E-17, -9.436970724121E-07 };
            int[] J0 = { 0, 1, -5, -4, -3, -2, -1, 2, 3 };
            double[] n0 = { -9.6927686500217, 10.086655968018, -0.005608791128302, 0.071452738081455, -0.40710498223928, 1.4240819171444, -4.383951131945, -0.28408632460772, 0.021268463753307 };
            double R = 0.461526; //kJ/(kg K)

            tau = 540 / T;
            g0_tau = 0;
            for (i = 0; i <= 8; i++)
            {
                g0_tau = g0_tau + n0[i] * J0[i] * Math.Pow(tau, J0[i] - 1);
            }
            gr_tau = 0;
            for (i = 0; i <= 42; i++)
            {
                gr_tau = gr_tau + nr[i] * Math.Pow(p, Ir[i]) * Jr[i] * Math.Pow(tau - 0.5, Jr[i] - 1);
            }
            return R * T * tau * (g0_tau + gr_tau);
        }

        public double u2_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //6 Equations for Region 2, Section. 6.1 Basic Equation
            //Table 11 and 12, Page 14 and 15
            int i = 0;
            double tau = 0;
            double g0_tau = 0;
            double g0_pi = 0;
            double gr_pi = 0;
            double gr_tau = 0;
            int[] Ir = { -1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24 };
            int[] Jr = { -0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58 };
            double[] nr = { -1.7731742473213E-03, -0.017834862292358, -0.045996013696365, -0.057581259083432, -0.05032527872793, -3.3032641670203E-05, -1.8948987516315E-04, -3.9392777243355E-03, -0.043797295650573, -2.6674547914087E-05, 2.0481737692309E-08, 4.3870667284435E-07, -3.227767723857E-05, -1.5033924542148E-03, -0.040668253562649, -7.8847309559367E-10, 1.2790717852285E-08, 4.8225372718507E-07, 2.2922076337661E-06, -1.6714766451061E-11, -2.1171472321355E-03, -23.895741934104, -5.905956432427E-18, -1.2621808899101E-06, -0.038946842435739, 1.1256211360459E-11, -8.2311340897998, 1.9809712802088E-08, 1.0406965210174E-19, -1.0234747095929E-13, -1.0018179379511E-09, -8.0882908646985E-11, 0.10693031879409, -0.33662250574171, 8.9185845355421E-25, 3.0629316876232E-13, -4.2002467698208E-06, -5.9056029685639E-26, 3.7826947613457E-06, -1.2768608934681E-15, 7.3087610595061E-29, 5.5414715350778E-17, -9.436970724121E-07 };
            int[] J0 = { 0, 1, -5, -4, -3, -2, -1, 2, 3 };
            double[] n0 = { -9.6927686500217, 10.086655968018, -0.005608791128302, 0.071452738081455, -0.40710498223928, 1.4240819171444, -4.383951131945, -0.28408632460772, 0.021268463753307 };
            double R = 0.461526; //kJ/(kg K)

            tau = 540 / T;
            g0_pi = 1 / p;
            g0_tau = 0;
            for (i = 0; i <= 8; i++)
            {
                g0_tau = g0_tau + n0[i] * J0[i] * Math.Pow(tau, J0[i] - 1);
            }
            gr_pi = 0;
            gr_tau = 0;
            for (i = 0; i <= 42; i++)
            {
                gr_pi = gr_pi + nr[i] * Ir[i] * Math.Pow(p, Ir[i] - 1) * Math.Pow(tau - 0.5, Jr[i]);
                gr_tau = gr_tau + nr[i] * Math.Pow(p, Ir[i]) * Jr[i] * Math.Pow(tau - 0.5, Jr[i] - 1);
            }
            return R * T * (tau * (g0_tau + gr_tau) - p * (g0_pi + gr_pi));
        }


        public double s2_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //6 Equations for Region 2, Section. 6.1 Basic Equation
            //Table 11 and 12, Page 14 and 15
            int i = 0;
            double tau = 0;
            double g0 = 0;
            double g0_tau = 0;
            double gr = 0;
            double gr_tau = 0;
            int[] Ir = { -1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24 };
            int[] Jr = { -0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58 };
            double[] nr = { -1.7731742473213E-03, -0.017834862292358, -0.045996013696365, -0.057581259083432, -0.05032527872793, -3.3032641670203E-05, -1.8948987516315E-04, -3.9392777243355E-03, -0.043797295650573, -2.6674547914087E-05, 2.0481737692309E-08, 4.3870667284435E-07, -3.227767723857E-05, -1.5033924542148E-03, -0.040668253562649, -7.8847309559367E-10, 1.2790717852285E-08, 4.8225372718507E-07, 2.2922076337661E-06, -1.6714766451061E-11, -2.1171472321355E-03, -23.895741934104, -5.905956432427E-18, -1.2621808899101E-06, -0.038946842435739, 1.1256211360459E-11, -8.2311340897998, 1.9809712802088E-08, 1.0406965210174E-19, -1.0234747095929E-13, -1.0018179379511E-09, -8.0882908646985E-11, 0.10693031879409, -0.33662250574171, 8.9185845355421E-25, 3.0629316876232E-13, -4.2002467698208E-06, -5.9056029685639E-26, 3.7826947613457E-06, -1.2768608934681E-15, 7.3087610595061E-29, 5.5414715350778E-17, -9.436970724121E-07 };
            int[] J0 = { 0, 1, -5, -4, -3, -2, -1, 2, 3 };
            double[] n0 = { -9.6927686500217, 10.086655968018, -0.005608791128302, 0.071452738081455, -0.40710498223928, 1.4240819171444, -4.383951131945, -0.28408632460772, 0.021268463753307 };
            double R = 0.461526; //kJ/(kg K)

            tau = 540 / T;
            g0 = Math.Log(p);
            g0_tau = 0;
            for (i = 0; i <= 8; i++)
            {
                g0 = g0 + n0[i] * Math.Pow(tau, J0[i]);
                g0_tau = g0_tau + n0[i] * J0[i] * Math.Pow(tau, J0[i] - 1);
            }
            gr = 0;
            gr_tau = 0;
            for (i = 0; i <= 42; i++)
            {
                gr = gr + nr[i] * Math.Pow(p, Ir[i]) * Math.Pow(tau - 0.5, Jr[i]);
                gr_tau = gr_tau + nr[i] * Math.Pow(p, Ir[i]) * Jr[i] * Math.Pow(tau - 0.5, Jr[i] - 1);
            }
            return R * (tau * (g0_tau + gr_tau) - (g0 + gr));
        }

        public double Cp2_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //6 Equations for Region 2, Section. 6.1 Basic Equation
            //Table 11 and 12, Page 14 and 15
            int i = 0;
            double tau = 0;
            double g0_tautau = 0;
            double gr_tautau = 0;
            int[] Ir = { -1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24 };
            int[] Jr = { -0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58 };
            double[] nr = { -1.7731742473213E-03, -0.017834862292358, -0.045996013696365, -0.057581259083432, -0.05032527872793, -3.3032641670203E-05, -1.8948987516315E-04, -3.9392777243355E-03, -0.043797295650573, -2.6674547914087E-05, 2.0481737692309E-08, 4.3870667284435E-07, -3.227767723857E-05, -1.5033924542148E-03, -0.040668253562649, -7.8847309559367E-10, 1.2790717852285E-08, 4.8225372718507E-07, 2.2922076337661E-06, -1.6714766451061E-11, -2.1171472321355E-03, -23.895741934104, -5.905956432427E-18, -1.2621808899101E-06, -0.038946842435739, 1.1256211360459E-11, -8.2311340897998, 1.9809712802088E-08, 1.0406965210174E-19, -1.0234747095929E-13, -1.0018179379511E-09, -8.0882908646985E-11, 0.10693031879409, -0.33662250574171, 8.9185845355421E-25, 3.0629316876232E-13, -4.2002467698208E-06, -5.9056029685639E-26, 3.7826947613457E-06, -1.2768608934681E-15, 7.3087610595061E-29, 5.5414715350778E-17, -9.436970724121E-07 };
            int[] J0 = { 0, 1, -5, -4, -3, -2, -1, 2, 3 };
            double[] n0 = { -9.6927686500217, 10.086655968018, -0.005608791128302, 0.071452738081455, -0.40710498223928, 1.4240819171444, -4.383951131945, -0.28408632460772, 0.021268463753307 };
            double R = 0.461526; //kJ/(kg K)

            tau = 540 / T;
            g0_tautau = 0;
            for (i = 0; i <= 8; i++)
            {
                g0_tautau = g0_tautau + n0[i] * J0[i] * (J0[i] - 1) * Math.Pow(tau, J0[i] - 2);
            }
            gr_tautau = 0;
            for (i = 0; i <= 42; i++)
            {
                gr_tautau = gr_tautau + nr[i] * Math.Pow(p, Ir[i]) * Jr[i] * (Jr[i] - 1) * Math.Pow(tau - 0.5, Jr[i] - 2);
            }
            return -R * Math.Pow(tau, 2) * (g0_tautau + gr_tautau);
        }

        public double Cv2_pT(double p, double T)
        {
            int i = 0;
            double tau = 0;
            double g0_tautau = 0;
            double gr_pi = 0;
            double gr_pitau = 0;
            double gr_pipi = 0;
            double gr_tautau = 0;
            int[] Ir = { -1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24 };
            int[] Jr = { -0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58 };
            double[] nr = { -1.7731742473213E-03, -0.017834862292358, -0.045996013696365, -0.057581259083432, -0.05032527872793, -3.3032641670203E-05, -1.8948987516315E-04, -3.9392777243355E-03, -0.043797295650573, -2.6674547914087E-05, 2.0481737692309E-08, 4.3870667284435E-07, -3.227767723857E-05, -1.5033924542148E-03, -0.040668253562649, -7.8847309559367E-10, 1.2790717852285E-08, 4.8225372718507E-07, 2.2922076337661E-06, -1.6714766451061E-11, -2.1171472321355E-03, -23.895741934104, -5.905956432427E-18, -1.2621808899101E-06, -0.038946842435739, 1.1256211360459E-11, -8.2311340897998, 1.9809712802088E-08, 1.0406965210174E-19, -1.0234747095929E-13, -1.0018179379511E-09, -8.0882908646985E-11, 0.10693031879409, -0.33662250574171, 8.9185845355421E-25, 3.0629316876232E-13, -4.2002467698208E-06, -5.9056029685639E-26, 3.7826947613457E-06, -1.2768608934681E-15, 7.3087610595061E-29, 5.5414715350778E-17, -9.436970724121E-07 };
            int[] J0 = { 0, 1, -5, -4, -3, -2, -1, 2, 3 };
            double[] n0 = { -9.6927686500217, 10.086655968018, -0.005608791128302, 0.071452738081455, -0.40710498223928, 1.4240819171444, -4.383951131945, -0.28408632460772, 0.021268463753307 };
            double R = 0.461526; //kJ/(kg K)

            tau = 540 / T;
            g0_tautau = 0;
            for (i = 0; i <= 8; i++)
            {
                g0_tautau = g0_tautau + n0[i] * J0[i] * (J0[i] - 1) * Math.Pow(tau, J0[i] - 2);
            }
            gr_pi = 0;
            gr_pitau = 0;
            gr_pipi = 0;
            gr_tautau = 0;
            for (i = 0; i <= 42; i++)
            {
                gr_pi = gr_pi + nr[i] * Ir[i] * Math.Pow(p, Ir[i] - 1) * Math.Pow(tau - 0.5, Jr[i]);
                gr_pipi = gr_pipi + nr[i] * Ir[i] * (Ir[i] - 1) * Math.Pow(p, Ir[i] - 2) * Math.Pow(tau - 0.5, Jr[i]);
                gr_pitau = gr_pitau + nr[i] * Ir[i] * Math.Pow(p, Ir[i] - 1) * Jr[i] * Math.Pow(tau - 0.5, Jr[i] - 1);
                gr_tautau = gr_tautau + nr[i] * Math.Pow(p, Ir[i]) * Jr[i] * (Jr[i] - 1) * Math.Pow(tau - 0.5, Jr[i] - 2);
            }
            return R * (-(Math.Pow(tau, 2) * (g0_tautau + gr_tautau)) - (Math.Pow(1 + p * gr_pi - tau * p * gr_pitau, 2)) / (1 - Math.Pow(p, 2) * gr_pipi));
        }

        public double w2_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //6 Equations for Region 2, Section. 6.1 Basic Equation
            //Table 11 and 12, Page 14 and 15
            int i = 0;
            double tau = 0;
            double g0_tautau = 0;
            double gr_pi = 0;
            double gr_pitau = 0;
            double gr_pipi = 0;
            double gr_tautau = 0;
            double R = 0.461526; //kJ/(kg K)
            int[] J0 = { 0, 1, -5, -4, -3, -2, -1, 2, 3 };
            double[] n0 = { -9.6927686500217, 10.086655968018, -0.005608791128302, 0.071452738081455, -0.40710498223928, 1.4240819171444, -4.383951131945, -0.28408632460772, 0.021268463753307 };
            int[] Ir = { -1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24 };
            int[] Jr = { -0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58 };
            double[] nr = { -1.7731742473213E-03, -0.017834862292358, -0.045996013696365, -0.057581259083432, -0.05032527872793, -3.3032641670203E-05, -1.8948987516315E-04, -3.9392777243355E-03, -0.043797295650573, -2.6674547914087E-05, 2.0481737692309E-08, 4.3870667284435E-07, -3.227767723857E-05, -1.5033924542148E-03, -0.040668253562649, -7.8847309559367E-10, 1.2790717852285E-08, 4.8225372718507E-07, 2.2922076337661E-06, -1.6714766451061E-11, -2.1171472321355E-03, -23.895741934104, -5.905956432427E-18, -1.2621808899101E-06, -0.038946842435739, 1.1256211360459E-11, -8.2311340897998, 1.9809712802088E-08, 1.0406965210174E-19, -1.0234747095929E-13, -1.0018179379511E-09, -8.0882908646985E-11, 0.10693031879409, -0.33662250574171, 8.9185845355421E-25, 3.0629316876232E-13, -4.2002467698208E-06, -5.9056029685639E-26, 3.7826947613457E-06, -1.2768608934681E-15, 7.3087610595061E-29, 5.5414715350778E-17, -9.436970724121E-07 };

            tau = 540 / T;
            g0_tautau = 0;
            for (i = 0; i <= 8; i++)
            {
                g0_tautau = g0_tautau + n0[i] * J0[i] * (J0[i] - 1) * Math.Pow(tau, J0[i] - 2);
            }
            gr_pi = 0;
            gr_pitau = 0;
            gr_pipi = 0;
            gr_tautau = 0;
            for (i = 0; i <= 42; i++)
            {
                gr_pi = gr_pi + nr[i] * Ir[i] * Math.Pow(p, Ir[i] - 1) * Math.Pow(tau - 0.5, Jr[i]);
                gr_pipi = gr_pipi + nr[i] * Ir[i] * (Ir[i] - 1) * Math.Pow(p, Ir[i] - 2) * Math.Pow(tau - 0.5, Jr[i]);
                gr_pitau = gr_pitau + nr[i] * Ir[i] * Math.Pow(p, Ir[i] - 1) * Jr[i] * Math.Pow(tau - 0.5, Jr[i] - 1);
                gr_tautau = gr_tautau + nr[i] * Math.Pow(p, Ir[i]) * Jr[i] * (Jr[i] - 1) * Math.Pow(tau - 0.5, Jr[i] - 2);
            }
            return Math.Pow(1000 * R * T * (1 + 2 * p * gr_pi + Math.Pow(p, 2) * Math.Pow(gr_pi, 2)) / ((1 - Math.Pow(p, 2) * gr_pipi) + Math.Pow(1 + p * gr_pi - tau * p * gr_pitau, 2) / (Math.Pow(tau, 2) * (g0_tautau + gr_tautau))), 0.5);
        }

        public double T2_ph(double p, double h)
        {
            double tempT2_ph = 0;
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //6 Equations for Region 2,6.3.1 The Backward Equations T( p, h ) for Subregions 2a, 2b, and 2c
            int sub_reg = 0;
            int i = 0;
            double Ts = 0;
            double hs = 0;


            if (p < 4)
            {
                sub_reg = 1;
            }
            else
            {
                if (p < (905.84278514723 - 0.67955786399241 * h + 1.2809002730136E-04 * Math.Pow(h, 2)))
                {
                    sub_reg = 2;
                }
                else
                {
                    sub_reg = 3;
                }
            }

            switch (sub_reg)
            {
                case 1:
                    //Subregion A
                    //Table 20, Eq 22, page 22
                    int[] Ji = { 0, 1, 2, 3, 7, 20, 0, 1, 2, 3, 7, 9, 11, 18, 44, 0, 2, 7, 36, 38, 40, 42, 44, 24, 44, 12, 32, 44, 32, 36, 42, 34, 44, 28 };
                    int[] Ii = { 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7 };
                    double[] ni = { 1089.8952318288, 849.51654495535, -107.81748091826, 33.153654801263, -7.4232016790248, 11.765048724356, 1.844574935579, -4.1792700549624, 6.2478196935812, -17.344563108114, -200.58176862096, 271.96065473796, -455.11318285818, 3091.9688604755, 252266.40357872, -6.1707422868339E-03, -0.31078046629583, 11.670873077107, 128127984.04046, -985549096.23276, 2822454697.3002, -3594897141.0703, 1722734991.3197, -13551.334240775, 12848734.66465, 1.3865724283226, 235988.32556514, -13105236.545054, 7399.9835474766, -551966.9703006, 3715408.5996233, 19127.72923966, -415351.64835634, -62.459855192507 };

                    Ts = 0;
                    hs = h / 2000;
                    for (i = 0; i <= 33; i++)
                    {
                        Ts = Ts + ni[i] * Math.Pow(p, Ii[i]) * Math.Pow(hs - 2.1, Ji[i]);
                    }
                    tempT2_ph = Ts;
                    break;
                case 2:
                    //Subregion B
                    //Table 21, Eq 23, page 23
                    Ji = new int[] { 0, 1, 2, 12, 18, 24, 28, 40, 0, 2, 6, 12, 18, 24, 28, 40, 2, 8, 18, 40, 1, 2, 12, 24, 2, 12, 18, 24, 28, 40, 18, 24, 40, 28, 2, 28, 1, 40 };
                    Ii = new int[] { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 7, 7, 9, 9 };
                    ni = new double[] { 1489.5041079516, 743.07798314034, -97.708318797837, 2.4742464705674, -0.63281320016026, 1.1385952129658, -0.47811863648625, 8.5208123431544E-03, 0.93747147377932, 3.3593118604916, 3.3809355601454, 0.16844539671904, 0.73875745236695, -0.47128737436186, 0.15020273139707, -0.002176411421975, -0.021810755324761, -0.10829784403677, -0.046333324635812, 7.1280351959551E-05, 1.1032831789999E-04, 1.8955248387902E-04, 3.0891541160537E-03, 1.3555504554949E-03, 2.8640237477456E-07, -1.0779857357512E-05, -7.6462712454814E-05, 1.4052392818316E-05, -3.1083814331434E-05, -1.0302738212103E-06, 2.821728163504E-07, 1.2704902271945E-06, 7.3803353468292E-08, -1.1030139238909E-08, -8.1456365207833E-14, -2.5180545682962E-11, -1.7565233969407E-18, 8.6934156344163E-15 };
                    Ts = 0;
                    hs = h / 2000;
                    for (i = 0; i <= 37; i++)
                    {
                        Ts = Ts + ni[i] * Math.Pow(p - 2, Ii[i]) * Math.Pow(hs - 2.6, Ji[i]);
                    }
                    tempT2_ph = Ts;
                    break;
                default:
                    //Subregion C
                    //Table 22, Eq 24, page 24



                    Ji = new int[] { 0, 4, 0, 2, 0, 2, 0, 1, 0, 2, 0, 1, 4, 8, 4, 0, 1, 4, 10, 12, 16, 20, 22 };
                    Ii = new int[] { -7, -7, -6, -6, -5, -5, -2, -2, -1, -1, 0, 0, 1, 1, 2, 6, 6, 6, 6, 6, 6, 6, 6 };
                    ni = new double[] { -3236839855524.2, 7326335090218.1, 358250899454.47, -583401318515.9, -10783068217.47, 20825544563.171, 610747.83564516, 859777.2253558, -25745.72360417, 31081.088422714, 1208.2315865936, 482.19755109255, 3.7966001272486, -10.842984880077, -0.04536417267666, 1.4559115658698E-13, 1.126159740723E-12, -1.7804982240686E-11, 1.2324579690832E-07, -1.1606921130984E-06, 2.7846367088554E-05, -5.9270038474176E-04, 1.2918582991878E-03 };
                    Ts = 0;
                    hs = h / 2000;
                    for (i = 0; i <= 22; i++)
                    {
                        Ts = Ts + ni[i] * Math.Pow(p + 25, Ii[i]) * Math.Pow(hs - 1.8, Ji[i]);
                    }
                    tempT2_ph = Ts;
                    break;
            }
            return tempT2_ph;
        }

        public double T2_ps(double p, double s)
        {
            double tempT2_ps = 0;
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //6 Equations for Region 2,6.3.2 The Backward Equations T( p, s ) for Subregions 2a, 2b, and 2c
            //Page 26
            int sub_reg = 0;
            int i = 0;
            double teta = 0;
            double sigma = 0;


            if (p < 4)
            {
                sub_reg = 1;
            }
            else
            {
                if (s < 5.85)
                {
                    sub_reg = 3;
                }
                else
                {
                    sub_reg = 2;
                }
            }
            switch (sub_reg)
            {
                case 1:
                    //Subregion A
                    //Table 25, Eq 25, page 26
                    double[] Ii = { -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.25, -1.25, -1.25, -1, -1, -1, -1, -1, -1, -0.75, -0.75, -0.5, -0.5, -0.5, -0.5, -0.25, -0.25, -0.25, -0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75, 0.75, 1, 1, 1.25, 1.25, 1.5, 1.5 };
                    int[] Ji = { -24, -23, -19, -13, -11, -10, -19, -15, -6, -26, -21, -17, -16, -9, -8, -15, -14, -26, -13, -9, -7, -27, -25, -11, -6, 1, 4, 8, 11, 0, 1, 5, 6, 10, 14, 16, 0, 4, 9, 17, 7, 18, 3, 15, 5, 18 };
                    double[] ni = { -392359.83861984, 515265.7382727, 40482.443161048, -321.93790923902, 96.961424218694, -22.867846371773, -449429.14124357, -5011.8336020166, 0.35684463560015, 44235.33584819, -13673.388811708, 421632.60207864, 22516.925837475, 474.42144865646, -149.31130797647, -197811.26320452, -23554.39947076, -19070.616302076, 55375.669883164, 3829.3691437363, -603.91860580567, 1936.3102620331, 4266.064369861, -5978.0638872718, -704.01463926862, 338.36784107553, 20.862786635187, 0.033834172656196, -4.3124428414893E-05, 166.53791356412, -139.86292055898, -0.78849547999872, 0.072132411753872, -5.9754839398283E-03, -1.2141358953904E-05, 2.3227096733871E-07, -10.538463566194, 2.0718925496502, -0.072193155260427, 2.074988708112E-07, -0.018340657911379, 2.9036272348696E-07, 0.21037527893619, 2.5681239729999E-04, -0.012799002933781, -8.2198102652018E-06 };


                    sigma = s / 2;
                    teta = 0;
                    for (i = 0; i <= 45; i++)
                    {
                        teta = teta + ni[i] * Math.Pow(p, Ii[i]) * Math.Pow(sigma - 2, Ji[i]);
                    }
                    tempT2_ps = teta;
                    break;
                case 2:
                    //Subregion B
                    //Table 26, Eq 26, page 27

                    Ii = new double[] { -6, -6, -5, -5, -4, -4, -4, -3, -3, -3, -3, -2, -2, -2, -2, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5 };
                    Ji = new int[] { 0, 11, 0, 11, 0, 1, 11, 0, 1, 11, 12, 0, 1, 6, 10, 0, 1, 5, 8, 9, 0, 1, 2, 4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 0, 1, 5, 0, 1, 3, 0, 1, 0, 1, 2 };
                    ni = new double[] { 316876.65083497, 20.864175881858, -398593.99803599, -21.816058518877, 223697.85194242, -2784.1703445817, 9.920743607148, -75197.512299157, 2970.8605951158, -3.4406878548526, 0.38815564249115, 17511.29508575, -1423.7112854449, 1.0943803364167, 0.89971619308495, -3375.9740098958, 471.62885818355, -1.9188241993679, 0.41078580492196, -0.33465378172097, 1387.0034777505, -406.63326195838, 41.72734715961, 2.1932549434532, -1.0320050009077, 0.35882943516703, 5.2511453726066E-03, 12.838916450705, -2.8642437219381, 0.56912683664855, -0.099962954584931, -3.2632037778459E-03, 2.3320922576723E-04, -0.1533480985745, 0.029072288239902, 3.7534702741167E-04, 1.7296691702411E-03, -3.8556050844504E-04, -3.5017712292608E-05, -1.4566393631492E-05, 5.6420857267269E-06, 4.1286150074605E-08, -2.0684671118824E-08, 1.6409393674725E-09 };


                    sigma = s / 0.7853;
                    teta = 0;
                    for (i = 0; i <= 43; i++)
                    {
                        teta = teta + ni[i] * Math.Pow(p, Ii[i]) * Math.Pow(10 - sigma, Ji[i]);
                    }
                    tempT2_ps = teta;
                    break;
                default:
                    //Subregion C
                    //Table 27, Eq 27, page 28

                    Ii = new double[] { -2, -2, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7, 7, 7 };
                    Ji = new int[] { 0, 1, 0, 0, 1, 2, 3, 0, 1, 3, 4, 0, 1, 2, 0, 1, 5, 0, 1, 4, 0, 1, 2, 0, 1, 0, 1, 3, 4, 5 };
                    ni = new double[] { 909.68501005365, 2404.566708842, -591.6232638713, 541.45404128074, -270.98308411192, 979.76525097926, -469.66772959435, 14.399274604723, -19.104204230429, 5.3299167111971, -21.252975375934, -0.3114733441376, 0.60334840894623, -0.042764839702509, 5.8185597255259E-03, -0.014597008284753, 5.6631175631027E-03, -7.6155864584577E-05, 2.2440342919332E-04, -1.2561095013413E-05, 6.3323132660934E-07, -2.0541989675375E-06, 3.6405370390082E-08, -2.9759897789215E-09, 1.0136618529763E-08, 5.9925719692351E-12, -2.0677870105164E-11, -2.0874278181886E-11, 1.0162166825089E-10, -1.6429828281347E-10 };


                    sigma = s / 2.9251;
                    teta = 0;
                    for (i = 0; i <= 29; i++)
                    {
                        teta = teta + ni[i] * Math.Pow(p, Ii[i]) * Math.Pow(2 - sigma, Ji[i]);
                    }
                    tempT2_ps = teta;
                    break;
            }
            return tempT2_ps;
        }

        public double p2_hs(double h, double s)
        {
            double tempp2_hs = 0;
            //Supplementary Release on Backward Equations for Pressure as a Function of Enthalpy and Entropy p(h,s) to the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam
            //Chapter 6:Backward Equations p(h,s) for Region 2
            int sub_reg = 0;
            int i = 0;
            double eta = 0;
          //  double teta = 0;
            double sigma = 0;
            double p = 0;


            if (h < (-3498.98083432139 + 2575.60716905876 * s - 421.073558227969 * Math.Pow(s, 2) + 27.6349063799944 * Math.Pow(s, 3)))
            {
                sub_reg = 1;
            }
            else
            {
                if (s < 5.85)
                {
                    sub_reg = 3;
                }
                else
                {
                    sub_reg = 2;
                }
            }
            switch (sub_reg)
            {
                case 1:
                    //Subregion A
                    //Table 6, Eq 3, page 8
                    int[] Ii = { 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 4, 5, 5, 6, 7 };
                    int[] Ji = { 1, 3, 6, 16, 20, 22, 0, 1, 2, 3, 5, 6, 10, 16, 20, 22, 3, 16, 20, 0, 2, 3, 6, 16, 16, 3, 16, 3, 1 };
                    double[] ni = { -1.82575361923032E-02, -0.125229548799536, 0.592290437320145, 6.04769706185122, 238.624965444474, -298.639090222922, 0.051225081304075, -0.437266515606486, 0.413336902999504, -5.16468254574773, -5.57014838445711, 12.8555037824478, 11.414410895329, -119.504225652714, -2847.7798596156, 4317.57846408006, 1.1289404080265, 1974.09186206319, 1516.12444706087, 1.41324451421235E-02, 0.585501282219601, -2.97258075863012, 5.94567314847319, -6236.56565798905, 9659.86235133332, 6.81500934948134, -6332.07286824489, -5.5891922446576, 4.00645798472063E-02 };


                    eta = h / 4200;
                    sigma = s / 12;
                    p = 0;
                    for (i = 0; i <= 28; i++)
                    {
                        p = p + ni[i] * Math.Pow(eta - 0.5, Ii[i]) * Math.Pow(sigma - 1.2, Ji[i]);
                    }
                    tempp2_hs = Math.Pow(p, 4) * 4;
                    break;
                case 2:
                    //Subregion B
                    //Table 7, Eq 4, page 9

                    Ii = new int[] { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6, 6, 7, 7, 8, 8, 8, 8, 12, 14 };
                    Ji = new int[] { 0, 1, 2, 4, 8, 0, 1, 2, 3, 5, 12, 1, 6, 18, 0, 1, 7, 12, 1, 16, 1, 12, 1, 8, 18, 1, 16, 1, 3, 14, 18, 10, 16 };
                    ni = new double[] { 8.01496989929495E-02, -0.543862807146111, 0.337455597421283, 8.9055545115745, 313.840736431485, 0.797367065977789, -1.2161697355624, 8.72803386937477, -16.9769781757602, -186.552827328416, 95115.9274344237, -18.9168510120494, -4334.0703719484, 543212633.012715, 0.144793408386013, 128.024559637516, -67230.9534071268, 33697238.0095287, -586.63419676272, -22140322476.9889, 1716.06668708389, -570817595.806302, -3121.09693178482, -2078413.8463301, 3056059461577.86, 3221.57004314333, 326810259797.295, -1441.04158934487, 410.694867802691, 109077066873.024, -24796465425889.3, 1888019068.65134, -123651009018773 };


                    eta = h / 4100;
                    sigma = s / 7.9;
                    for (i = 0; i <= 32; i++)
                    {
                        p = p + ni[i] * Math.Pow(eta - 0.6, Ii[i]) * Math.Pow(sigma - 1.01, Ji[i]);
                    }
                    tempp2_hs = Math.Pow(p, 4) * 100;
                    break;
                default:
                    //Subregion C
                    //Table 8, Eq 5, page 10


                    Ii = new int[] { 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5, 6, 6, 10, 12, 16 };
                    Ji = new int[] { 0, 1, 2, 3, 4, 8, 0, 2, 5, 8, 14, 2, 3, 7, 10, 18, 0, 5, 8, 16, 18, 18, 1, 4, 6, 14, 8, 18, 7, 7, 10 };
                    ni = new double[] { 0.112225607199012, -3.39005953606712, -32.0503911730094, -197.5973051049, -407.693861553446, 13294.3775222331, 1.70846839774007, 37.3694198142245, 3581.44365815434, 423014.446424664, -751071025.760063, 52.3446127607898, -228.351290812417, -960652.417056937, -80705929.2526074, 1626980172256.69, 0.772465073604171, 46392.9973837746, -13731788.5134128, 1704703926305.12, -25110462818730.8, 31774883083552, 53.8685623675312, -55308.9094625169, -1028615.22421405, 2042494187562.34, 273918446.626977, -2.63963146312685E+15, -1078908541.08088, -29649262098.0124, -1.11754907323424E+15 };

                    eta = h / 3500;
                    sigma = s / 5.9;
                    for (i = 0; i <= 30; i++)
                    {
                        p = p + ni[i] * Math.Pow(eta - 0.7, Ii[i]) * Math.Pow(sigma - 1.1, Ji[i]);
                    }
                    tempp2_hs = Math.Pow(p, 4) * 100;
                    break;
            }
            return tempp2_hs;
        }

        public double T2_prho(double p, double rho)
        {
            //Solve by iteration. Observe that fo low temperatures this equation has 2 solutions.
            //Solve with half interval method
            double Low_Bound = 0;
            double High_Bound = 0;
            double rhos = 0;
            double Ts = 0;

            if (p < 16.5292)
            {
                Low_Bound = T4_p(p);
            }
            else
            {
                Low_Bound = B23T_p(p);
            }
            High_Bound = 1073.15;
            while (Math.Abs(rho - rhos) > 0.000001)
            {
                Ts = (Low_Bound + High_Bound) / 2;
                rhos = 1 / v2_pT(p, Ts);
                if (rhos < rho)
                {
                    High_Bound = Ts;
                }
                else
                {
                    Low_Bound = Ts;
                }
            }
            return Ts;
        }


        //***********************************************************************************************************
        //*2.3 Functions for region 3

        public double p3_rhoT(double rho, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //7 Basic Equation for Region 3, Section. 6.1 Basic Equation
            //Table 30 and 31, Page 30 and 31
            int i = 0;

            double delta = 0;
            double tau = 0;
            double fidelta = 0;
            double R = 0.461526;
            double tc = 647.096;
       //     double pc = 22.064;
            double rhoc = 322;

            int[] Ii = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11 };
            int[] Ji = { 0, 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26 };
            double[] ni = { 1.0658070028513, -15.732845290239, 20.944396974307, -7.6867707878716, 2.6185947787954, -2.808078114862, 1.2053369696517, -8.4566812812502E-03, -1.2654315477714, -1.1524407806681, 0.88521043984318, -0.64207765181607, 0.38493460186671, -0.85214708824206, 4.8972281541877, -3.0502617256965, 0.039420536879154, 0.12558408424308, -0.2799932969871, 1.389979956946, -2.018991502357, -8.2147637173963E-03, -0.47596035734923, 0.0439840744735, -0.44476435428739, 0.90572070719733, 0.70522450087967, 0.10770512626332, -0.32913623258954, -0.50871062041158, -0.022175400873096, 0.094260751665092, 0.16436278447961, -0.013503372241348, -0.014834345352472, 5.7922953628084E-04, 3.2308904703711E-03, 8.0964802996215E-05, -1.6557679795037E-04, -4.4923899061815E-05 };


            delta = rho / rhoc;
            tau = tc / T;
            fidelta = 0;
            for (i = 1; i <= 39; i++)
            {
                fidelta = fidelta + ni[i] * Ii[i] * Math.Pow(delta, Ii[i] - 1) * Math.Pow(tau, Ji[i]);
            }
            fidelta = fidelta + ni[0] / delta;
            return rho * R * T * delta * fidelta / 1000;
        }

        public double u3_rhoT(double rho, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //7 Basic Equation for Region 3, Section. 6.1 Basic Equation
            //Table 30 and 31, Page 30 and 31
            int i = 0;

            double delta = 0;
            double tau = 0;
            double fitau = 0;
            double R = 0.461526;
            double tc = 647.096;
         //   double pc = 22.064;
            double rhoc = 322;
            int[] Ii = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11 };
            int[] Ji = { 0, 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26 };
            double[] ni = { 1.0658070028513, -15.732845290239, 20.944396974307, -7.6867707878716, 2.6185947787954, -2.808078114862, 1.2053369696517, -8.4566812812502E-03, -1.2654315477714, -1.1524407806681, 0.88521043984318, -0.64207765181607, 0.38493460186671, -0.85214708824206, 4.8972281541877, -3.0502617256965, 0.039420536879154, 0.12558408424308, -0.2799932969871, 1.389979956946, -2.018991502357, -8.2147637173963E-03, -0.47596035734923, 0.0439840744735, -0.44476435428739, 0.90572070719733, 0.70522450087967, 0.10770512626332, -0.32913623258954, -0.50871062041158, -0.022175400873096, 0.094260751665092, 0.16436278447961, -0.013503372241348, -0.014834345352472, 5.7922953628084E-04, 3.2308904703711E-03, 8.0964802996215E-05, -1.6557679795037E-04, -4.4923899061815E-05 };

            delta = rho / rhoc;
            tau = tc / T;
            fitau = 0;
            for (i = 1; i <= 39; i++)
            {
                fitau = fitau + ni[i] * Math.Pow(delta, Ii[i]) * Ji[i] * Math.Pow(tau, Ji[i] - 1);
            }
            return R * T * (tau * fitau);
        }

        public double h3_rhoT(double rho, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //7 Basic Equation for Region 3, Section. 6.1 Basic Equation
            //Table 30 and 31, Page 30 and 31
            int i = 0;
            double delta = 0;
            double tau = 0;
            double fidelta = 0;
            double fitau = 0;
            double R = 0.461526;
            double tc = 647.096;
         //   double pc = 22.064;
            double rhoc = 322;

            int[] Ii = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11 };
            int[] Ji = { 0, 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26 };
            double[] ni = { 1.0658070028513, -15.732845290239, 20.944396974307, -7.6867707878716, 2.6185947787954, -2.808078114862, 1.2053369696517, -8.4566812812502E-03, -1.2654315477714, -1.1524407806681, 0.88521043984318, -0.64207765181607, 0.38493460186671, -0.85214708824206, 4.8972281541877, -3.0502617256965, 0.039420536879154, 0.12558408424308, -0.2799932969871, 1.389979956946, -2.018991502357, -8.2147637173963E-03, -0.47596035734923, 0.0439840744735, -0.44476435428739, 0.90572070719733, 0.70522450087967, 0.10770512626332, -0.32913623258954, -0.50871062041158, -0.022175400873096, 0.094260751665092, 0.16436278447961, -0.013503372241348, -0.014834345352472, 5.7922953628084E-04, 3.2308904703711E-03, 8.0964802996215E-05, -1.6557679795037E-04, -4.4923899061815E-05 };

            delta = rho / rhoc;
            tau = tc / T;
            fidelta = 0;
            fitau = 0;
            for (i = 1; i <= 39; i++)
            {
                fidelta = fidelta + ni[i] * Ii[i] * Math.Pow(delta, Ii[i] - 1) * Math.Pow(tau, Ji[i]);
                fitau = fitau + ni[i] * Math.Pow(delta, Ii[i]) * Ji[i] * Math.Pow(tau, Ji[i] - 1);
            }
            fidelta = fidelta + ni[0] / delta;
            return R * T * (tau * fitau + delta * fidelta);
        }

        public double s3_rhoT(double rho, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //7 Basic Equation for Region 3, Section. 6.1 Basic Equation
            //Table 30 and 31, Page 30 and 31
            int i = 0;
            double fi = 0;
            double delta = 0;
            double tau = 0;
            double fitau = 0;
            double R = 0.461526;
            double tc = 647.096;
         //   double pc = 22.064;
            double rhoc = 322;
            int[] Ii = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11 };
            int[] Ji = { 0, 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26 };
            double[] ni = { 1.0658070028513, -15.732845290239, 20.944396974307, -7.6867707878716, 2.6185947787954, -2.808078114862, 1.2053369696517, -8.4566812812502E-03, -1.2654315477714, -1.1524407806681, 0.88521043984318, -0.64207765181607, 0.38493460186671, -0.85214708824206, 4.8972281541877, -3.0502617256965, 0.039420536879154, 0.12558408424308, -0.2799932969871, 1.389979956946, -2.018991502357, -8.2147637173963E-03, -0.47596035734923, 0.0439840744735, -0.44476435428739, 0.90572070719733, 0.70522450087967, 0.10770512626332, -0.32913623258954, -0.50871062041158, -0.022175400873096, 0.094260751665092, 0.16436278447961, -0.013503372241348, -0.014834345352472, 5.7922953628084E-04, 3.2308904703711E-03, 8.0964802996215E-05, -1.6557679795037E-04, -4.4923899061815E-05 };

            delta = rho / rhoc;
            tau = tc / T;
            fi = 0;
            fitau = 0;
            for (i = 1; i <= 39; i++)
            {
                fi = fi + ni[i] * Math.Pow(delta, Ii[i]) * Math.Pow(tau, Ji[i]);
                fitau = fitau + ni[i] * Math.Pow(delta, Ii[i]) * Ji[i] * Math.Pow(tau, Ji[i] - 1);
            }
            fi = fi + ni[0] * Math.Log(delta);
            return R * (tau * fitau - fi);
        }

        public double Cp3_rhoT(double rho, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //7 Basic Equation for Region 3, Section. 6.1 Basic Equation
            //Table 30 and 31, Page 30 and 31
            int i = 0;

            double fideltatau = 0;
         //   double fi = 0;
            double delta = 0;
            double tau = 0;
            double fitautau = 0;
            double fidelta = 0;
           // double fideltatautau = 0;
            double fideltadelta = 0;
            double R = 0.461526;
            double tc = 647.096;
          //  double pc = 22.064;
            double rhoc = 322;

            int[] Ii = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11 };
            int[] Ji = { 0, 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26 };
            double[] ni = { 1.0658070028513, -15.732845290239, 20.944396974307, -7.6867707878716, 2.6185947787954, -2.808078114862, 1.2053369696517, -8.4566812812502E-03, -1.2654315477714, -1.1524407806681, 0.88521043984318, -0.64207765181607, 0.38493460186671, -0.85214708824206, 4.8972281541877, -3.0502617256965, 0.039420536879154, 0.12558408424308, -0.2799932969871, 1.389979956946, -2.018991502357, -8.2147637173963E-03, -0.47596035734923, 0.0439840744735, -0.44476435428739, 0.90572070719733, 0.70522450087967, 0.10770512626332, -0.32913623258954, -0.50871062041158, -0.022175400873096, 0.094260751665092, 0.16436278447961, -0.013503372241348, -0.014834345352472, 5.7922953628084E-04, 3.2308904703711E-03, 8.0964802996215E-05, -1.6557679795037E-04, -4.4923899061815E-05 };


            delta = rho / rhoc;
            tau = tc / T;
            fitautau = 0;
            fidelta = 0;
            fideltatau = 0;
            fideltadelta = 0;
            for (i = 1; i <= 39; i++)
            {
                fitautau = fitautau + ni[i] * Math.Pow(delta, Ii[i]) * Ji[i] * (Ji[i] - 1) * Math.Pow(tau, Ji[i] - 2);
                fidelta = fidelta + ni[i] * Ii[i] * Math.Pow(delta, Ii[i] - 1) * Math.Pow(tau, Ji[i]);
                fideltatau = fideltatau + ni[i] * Ii[i] * Math.Pow(delta, Ii[i] - 1) * Ji[i] * Math.Pow(tau, Ji[i] - 1);
                fideltadelta = fideltadelta + ni[i] * Ii[i] * (Ii[i] - 1) * Math.Pow(delta, Ii[i] - 2) * Math.Pow(tau, Ji[i]);
            }
            fidelta = fidelta + ni[0] / delta;
            fideltadelta = fideltadelta - ni[0] / (Math.Pow(delta, 2));
            return R * (-(Math.Pow(tau, 2) * fitautau) + Math.Pow(delta * fidelta - delta * tau * fideltatau, 2) / (2 * delta * fidelta + Math.Pow(delta, 2) * fideltadelta));
        }

        public double Cv3_rhoT(double rho, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //7 Basic Equation for Region 3, Section. 6.1 Basic Equation
            //Table 30 and 31, Page 30 and 31
            int i = 0;

         //   double fi = 0;
            double delta = 0;
            double tau = 0;
            double fitautau = 0;
            double R = 0.461526;
            double tc = 647.096;
          //  double pc = 22.064;
            double rhoc = 322;

            int[] Ii = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11 };
            int[] Ji = { 0, 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26 };
            double[] ni = { 1.0658070028513, -15.732845290239, 20.944396974307, -7.6867707878716, 2.6185947787954, -2.808078114862, 1.2053369696517, -8.4566812812502E-03, -1.2654315477714, -1.1524407806681, 0.88521043984318, -0.64207765181607, 0.38493460186671, -0.85214708824206, 4.8972281541877, -3.0502617256965, 0.039420536879154, 0.12558408424308, -0.2799932969871, 1.389979956946, -2.018991502357, -8.2147637173963E-03, -0.47596035734923, 0.0439840744735, -0.44476435428739, 0.90572070719733, 0.70522450087967, 0.10770512626332, -0.32913623258954, -0.50871062041158, -0.022175400873096, 0.094260751665092, 0.16436278447961, -0.013503372241348, -0.014834345352472, 5.7922953628084E-04, 3.2308904703711E-03, 8.0964802996215E-05, -1.6557679795037E-04, -4.4923899061815E-05 };


            delta = rho / rhoc;
            tau = tc / T;
            fitautau = 0;
            for (i = 1; i <= 39; i++)
            {
                fitautau = fitautau + ni[i] * Math.Pow(delta, Ii[i]) * Ji[i] * (Ji[i] - 1) * Math.Pow(tau, Ji[i] - 2);
            }
            return R * -(tau * tau * fitautau);
        }

        public double w3_rhoT(double rho, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //7 Basic Equation for Region 3, Section. 6.1 Basic Equation
            //Table 30 and 31, Page 30 and 31
            int i = 0;

          //  double fi = 0;
            double delta = 0;
            double tau = 0;
            double fitautau = 0;
            double fidelta = 0;
            double fideltatau = 0;
            double fideltadelta = 0;
            double R = 0.461526;
            double tc = 647.096;
         //   double pc = 22.064;
            double rhoc = 322;
            int[] Ii = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11 };
            int[] Ji = { 0, 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26 };
            double[] ni = { 1.0658070028513, -15.732845290239, 20.944396974307, -7.6867707878716, 2.6185947787954, -2.808078114862, 1.2053369696517, -8.4566812812502E-03, -1.2654315477714, -1.1524407806681, 0.88521043984318, -0.64207765181607, 0.38493460186671, -0.85214708824206, 4.8972281541877, -3.0502617256965, 0.039420536879154, 0.12558408424308, -0.2799932969871, 1.389979956946, -2.018991502357, -8.2147637173963E-03, -0.47596035734923, 0.0439840744735, -0.44476435428739, 0.90572070719733, 0.70522450087967, 0.10770512626332, -0.32913623258954, -0.50871062041158, -0.022175400873096, 0.094260751665092, 0.16436278447961, -0.013503372241348, -0.014834345352472, 5.7922953628084E-04, 3.2308904703711E-03, 8.0964802996215E-05, -1.6557679795037E-04, -4.4923899061815E-05 };
            delta = rho / rhoc;
            tau = tc / T;
            fitautau = 0;
            fidelta = 0;
            fideltatau = 0;
            fideltadelta = 0;
            for (i = 1; i <= 39; i++)
            {
                fitautau = fitautau + ni[i] * Math.Pow(delta, Ii[i]) * Ji[i] * (Ji[i] - 1) * Math.Pow(tau, Ji[i] - 2);
                fidelta = fidelta + ni[i] * Ii[i] * Math.Pow(delta, Ii[i] - 1) * Math.Pow(tau, Ji[i]);
                fideltatau = fideltatau + ni[i] * Ii[i] * Math.Pow(delta, Ii[i] - 1) * Ji[i] * Math.Pow(tau, Ji[i] - 1);
                fideltadelta = fideltadelta + ni[i] * Ii[i] * (Ii[i] - 1) * Math.Pow(delta, Ii[i] - 2) * Math.Pow(tau, Ji[i]);
            }
            fidelta = fidelta + ni[0] / delta;
            fideltadelta = fideltadelta - ni[0] / (Math.Pow(delta, 2));
            return Math.Pow(1000 * R * T * (2 * delta * fidelta + Math.Pow(delta, 2) * fideltadelta - Math.Pow(delta * fidelta - delta * tau * fideltatau, 2) / (Math.Pow(tau, 2) * fitautau)), 0.5);
        }


        public double T3_ph(double p, double h)
        {
            double tempT3_ph = 0;
            //Revised Supplementary Release on Backward Equations for the Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam
            //2004
            //Section 3.3 Backward Equations T(p,h) and v(p,h) for Subregions 3a and 3b
            //Boundary equation, Eq 1 Page 5
            int i = 0;

            double h3ab = 0;
            double ps = 0;
            double hs = 0;
            double Ts = 0;
         //   double R = 0.461526;
         //   double tc = 647.096;
         //   double pc = 22.064;
        //    double rhoc = 322;
            h3ab = (2014.64004206875 + 3.74696550136983 * p - 2.19921901054187E-02 * Math.Pow(p, 2) + 8.7513168600995E-05 * Math.Pow(p, 3));
            if (h < h3ab)
            {
                //Subregion 3a
                //Eq 2, Table 3, Page 7

                int[] Ii = { -12, -12, -12, -12, -12, -12, -12, -12, -10, -10, -10, -8, -8, -8, -8, -5, -3, -2, -2, -2, -1, -1, 0, 0, 1, 3, 3, 4, 4, 10, 12 };
                int[] Ji = { 0, 1, 2, 6, 14, 16, 20, 22, 1, 5, 12, 0, 2, 4, 10, 2, 0, 1, 3, 4, 0, 2, 0, 1, 1, 0, 1, 0, 3, 4, 5 };
                double[] ni = { -1.33645667811215E-07, 4.55912656802978E-06, -1.46294640700979E-05, 6.3934131297008E-03, 372.783927268847, -7186.54377460447, 573494.7521034, -2675693.29111439, -3.34066283302614E-05, -2.45479214069597E-02, 47.8087847764996, 7.64664131818904E-06, 1.28350627676972E-03, 1.71219081377331E-02, -8.51007304583213, -1.36513461629781E-02, -3.84460997596657E-06, 3.37423807911655E-03, -0.551624873066791, 0.72920227710747, -9.92522757376041E-03, -0.119308831407288, 0.793929190615421, 0.454270731799386, 0.20999859125991, -6.42109823904738E-03, -0.023515586860454, 2.52233108341612E-03, -7.64885133368119E-03, 1.36176427574291E-02, -1.33027883575669E-02 };

                ps = p / 100;
                hs = h / 2300;
                Ts = 0;
                for (i = 0; i <= 30; i++)
                {
                    Ts = Ts + ni[i] * Math.Pow(ps + 0.24, Ii[i]) * Math.Pow(hs - 0.615, Ji[i]);
                }
                tempT3_ph = Ts * 760;
            }
            else
            {
                //Subregion 3b
                //Eq 3, Table 4, Page 7,8
                int[] Ii = { -12, -12, -10, -10, -10, -10, -10, -8, -8, -8, -8, -8, -6, -6, -6, -4, -4, -3, -2, -2, -1, -1, -1, -1, -1, -1, 0, 0, 1, 3, 5, 6, 8 };
                int[] Ji = { 0, 1, 0, 1, 5, 10, 12, 0, 1, 2, 4, 10, 0, 1, 2, 0, 1, 5, 0, 4, 2, 4, 6, 10, 14, 16, 0, 2, 1, 1, 1, 1, 1 };
                double[] ni = { 3.2325457364492E-05, -1.27575556587181E-04, -4.75851877356068E-04, 1.56183014181602E-03, 0.105724860113781, -85.8514221132534, 724.140095480911, 2.96475810273257E-03, -5.92721983365988E-03, -1.26305422818666E-02, -0.115716196364853, 84.9000969739595, -1.08602260086615E-02, 1.54304475328851E-02, 7.50455441524466E-02, 2.52520973612982E-02, -6.02507901232996E-02, -3.07622221350501, -5.74011959864879E-02, 5.03471360939849, -0.925081888584834, 3.91733882917546, -77.314600713019, 9493.08762098587, -1410437.19679409, 8491662.30819026, 0.861095729446704, 0.32334644281172, 0.873281936020439, -0.436653048526683, 0.286596714529479, -0.131778331276228, 6.76682064330275E-03 };



                hs = h / 2800;
                ps = p / 100;
                Ts = 0;
                for (i = 0; i <= 32; i++)
                {
                    Ts = Ts + ni[i] * Math.Pow(ps + 0.298, Ii[i]) * Math.Pow(hs - 0.72, Ji[i]);
                }
                tempT3_ph = Ts * 860;
            }
            return tempT3_ph;
        }

        public double v3_ph(double p, double h)
        {
            double tempv3_ph = 0;
            //Revised Supplementary Release on Backward Equations for the Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam
            //2004
            //Section 3.3 Backward Equations T(p,h) and v(p,h) for Subregions 3a and 3b
            //Boundary equation, Eq 1 Page 5
            int i = 0;
            double h3ab = 0;
            double ps = 0;
            double hs = 0;
            double vs = 0;
      //      double R = 0.461526;
        //    double tc = 647.096;
         //   double pc = 22.064;
       //     double rhoc = 322;
            h3ab = (2014.64004206875 + 3.74696550136983 * p - 2.19921901054187E-02 * Math.Pow(p, 2) + 8.7513168600995E-05 * Math.Pow(p, 3));
            if (h < h3ab)
            {
                //Subregion 3a
                //Eq 4, Table 6, Page 9
                int[] Ii = { -12, -12, -12, -12, -10, -10, -10, -8, -8, -6, -6, -6, -4, -4, -3, -2, -2, -1, -1, -1, -1, 0, 0, 1, 1, 1, 2, 2, 3, 4, 5, 8 };
                int[] Ji = { 6, 8, 12, 18, 4, 7, 10, 5, 12, 3, 4, 22, 2, 3, 7, 3, 16, 0, 1, 2, 3, 0, 1, 0, 1, 2, 0, 2, 0, 2, 2, 2 };
                double[] ni = { 5.29944062966028E-03, -0.170099690234461, 11.1323814312927, -2178.98123145125, -5.06061827980875E-04, 0.556495239685324, -9.43672726094016, -0.297856807561527, 93.9353943717186, 1.92944939465981E-02, 0.421740664704763, -3689141.2628233, -7.37566847600639E-03, -0.354753242424366, -1.99768169338727, 1.15456297059049, 5683.6687581596, 8.08169540124668E-03, 0.172416341519307, 1.04270175292927, -0.297691372792847, 0.560394465163593, 0.275234661176914, -0.148347894866012, -6.51142513478515E-02, -2.92468715386302, 6.64876096952665E-02, 3.52335014263844, -1.46340792313332E-02, -2.24503486668184, 1.10533464706142, -4.08757344495612E-02 };


                ps = p / 100;
                hs = h / 2100;
                vs = 0;
                for (i = 0; i <= 31; i++)
                {
                    vs = vs + ni[i] * Math.Pow(ps + 0.128, Ii[i]) * Math.Pow(hs - 0.727, Ji[i]);
                }
                tempv3_ph = vs * 0.0028;
            }
            else
            {
                //Subregion 3b
                //Eq 5, Table 7, Page 9
                int[] Ii = { -12, -12, -8, -8, -8, -8, -8, -8, -6, -6, -6, -6, -6, -6, -4, -4, -4, -3, -3, -2, -2, -1, -1, -1, -1, 0, 1, 1, 2, 2 };
                int[] Ji = { 0, 1, 0, 1, 3, 6, 7, 8, 0, 1, 2, 5, 6, 10, 3, 6, 10, 0, 2, 1, 2, 0, 1, 4, 5, 0, 0, 1, 2, 6 };
                double[] ni = { -2.25196934336318E-09, 1.40674363313486E-08, 2.3378408528056E-06, -3.31833715229001E-05, 1.07956778514318E-03, -0.271382067378863, 1.07202262490333, -0.853821329075382, -2.15214194340526E-05, 7.6965608822273E-04, -4.31136580433864E-03, 0.453342167309331, -0.507749535873652, -100.475154528389, -0.219201924648793, -3.21087965668917, 607.567815637771, 5.57686450685932E-04, 0.18749904002955, 9.05368030448107E-03, 0.285417173048685, 3.29924030996098E-02, 0.239897419685483, 4.82754995951394, -11.8035753702231, 0.169490044091791, -1.79967222507787E-02, 3.71810116332674E-02, -5.36288335065096E-02, 1.6069710109252 };
                ps = p / 100;
                hs = h / 2800;
                vs = 0;
                for (i = 0; i <= 29; i++)
                {
                    vs = vs + ni[i] * Math.Pow(ps + 0.0661, Ii[i]) * Math.Pow(hs - 0.72, Ji[i]);
                }
                tempv3_ph = vs * 0.0088;
            }
            return tempv3_ph;
        }


        public double T3_ps(double p, double s)
        {
            double tempT3_ps = 0;
            //Revised Supplementary Release on Backward Equations for the Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam
            //2004
            //3.4 Backward Equations T(p,s) and v(p,s) for Subregions 3a and 3b
            //Boundary equation, Eq 6 Page 11
            int i = 0;

            double ps = 0;
            double sigma = 0;
            double teta = 0;
       //     double R = 0.461526;
       //     double tc = 647.096;
        //    double pc = 22.064;
        //    double rhoc = 322;

            if (s <= 4.41202148223476)
            {
                //Subregion 3a
                //Eq 6, Table 10, Page 11

                int[] Ii = { -12, -12, -10, -10, -10, -10, -8, -8, -8, -8, -6, -6, -6, -5, -5, -5, -4, -4, -4, -2, -2, -1, -1, 0, 0, 0, 1, 2, 2, 3, 8, 8, 10 };
                int[] Ji = { 28, 32, 4, 10, 12, 14, 5, 7, 8, 28, 2, 6, 32, 0, 14, 32, 6, 10, 36, 1, 4, 1, 6, 0, 1, 4, 0, 0, 3, 2, 0, 1, 2 };
                double[] ni = { 1500420082.63875, -159397258480.424, 5.02181140217975E-04, -67.2057767855466, 1450.58545404456, -8238.8953488889, -0.154852214233853, 11.2305046746695, -29.7000213482822, 43856513263.5495, 1.37837838635464E-03, -2.97478527157462, 9717779473494.13, -5.71527767052398E-05, 28830.794977842, -74442828926270.3, 12.8017324848921, -368.275545889071, 6.64768904779177E+15, 0.044935925195888, -4.22897836099655, -0.240614376434179, -4.74341365254924, 0.72409399912611, 0.923874349695897, 3.99043655281015, 3.84066651868009E-02, -3.59344365571848E-03, -0.735196448821653, 0.188367048396131, 1.41064266818704E-04, -2.57418501496337E-03, 1.23220024851555E-03 };

                sigma = s / 4.4;
                ps = p / 100;
                teta = 0;
                for (i = 0; i <= 32; i++)
                {
                    teta = teta + ni[i] * Math.Pow(ps + 0.24, Ii[i]) * Math.Pow(sigma - 0.703, Ji[i]);
                }
                tempT3_ps = teta * 760;
            }
            else
            {
                //Subregion 3b
                //Eq 7, Table 11, Page 11
                int[] Ii = { -12, -12, -12, -12, -8, -8, -8, -6, -6, -6, -5, -5, -5, -5, -5, -4, -3, -3, -2, 0, 2, 3, 4, 5, 6, 8, 12, 14 };
                int[] Ji = { 1, 3, 4, 7, 0, 1, 3, 0, 2, 4, 0, 1, 2, 4, 6, 12, 1, 6, 2, 0, 1, 1, 0, 24, 0, 3, 1, 2 };
                double[] ni = { 0.52711170160166, -40.1317830052742, 153.020073134484, -2247.99398218827, -0.193993484669048, -1.40467557893768, 42.6799878114024, 0.752810643416743, 22.6657238616417, -622.873556909932, -0.660823667935396, 0.841267087271658, -25.3717501764397, 485.708963532948, 880.531517490555, 2650155.92794626, -0.359287150025783, -656.991567673753, 2.41768149185367, 0.856873461222588, 0.655143675313458, -0.213535213206406, 5.62974957606348E-03, -316955725450471, -6.99997000152457E-04, 1.19845803210767E-02, 1.93848122022095E-05, -2.15095749182309E-05 };

                sigma = s / 5.3;
                ps = p / 100;
                teta = 0;
                for (i = 0; i <= 27; i++)
                {
                    teta = teta + ni[i] * Math.Pow(ps + 0.76, Ii[i]) * Math.Pow(sigma - 0.818, Ji[i]);
                }
                tempT3_ps = teta * 860;
            }
            return tempT3_ps;
        }

        public double v3_ps(double p, double s)
        {
            double tempv3_ps = 0;
            //Revised Supplementary Release on Backward Equations for the Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam
            //2004
            //3.4 Backward Equations T(p,s) and v(p,s) for Subregions 3a and 3b
            //Boundary equation, Eq 6 Page 11
            int i = 0;

            double ps = 0;
            double sigma = 0;
            double omega = 0;
           // double R = 0.461526;
         //   double tc = 647.096;
        //    double pc = 22.064;
        //    double rhoc = 322;
            if (s <= 4.41202148223476)
            {
                //Subregion 3a
                //Eq 8, Table 13, Page 14
                int[] Ii = { -12, -12, -12, -10, -10, -10, -10, -8, -8, -8, -8, -6, -5, -4, -3, -3, -2, -2, -1, -1, 0, 0, 0, 1, 2, 4, 5, 6 };
                int[] Ji = { 10, 12, 14, 4, 8, 10, 20, 5, 6, 14, 16, 28, 1, 5, 2, 4, 3, 8, 1, 2, 0, 1, 3, 0, 0, 2, 2, 0 };
                double[] ni = { 79.5544074093975, -2382.6124298459, 17681.3100617787, -1.10524727080379E-03, -15.3213833655326, 297.544599376982, -35031520.6871242, 0.277513761062119, -0.523964271036888, -148011.182995403, 1600148.99374266, 1708023226634.27, 2.46866996006494E-04, 1.6532608479798, -0.118008384666987, 2.537986423559, 0.965127704669424, -28.2172420532826, 0.203224612353823, 1.10648186063513, 0.52612794845128, 0.277000018736321, 1.08153340501132, -7.44127885357893E-02, 1.64094443541384E-02, -6.80468275301065E-02, 0.025798857610164, -1.45749861944416E-04 };
                ps = p / 100;
                sigma = s / 4.4;
                omega = 0;
                for (i = 0; i <= 27; i++)
                {
                    omega = omega + ni[i] * Math.Pow(ps + 0.187, Ii[i]) * Math.Pow(sigma - 0.755, Ji[i]);
                }
                tempv3_ps = omega * 0.0028;
            }
            else
            {
                //Subregion 3b
                //Eq 9, Table 14, Page 14
                int[] Ii = { -12, -12, -12, -12, -12, -12, -10, -10, -10, -10, -8, -5, -5, -5, -4, -4, -4, -4, -3, -2, -2, -2, -2, -2, -2, 0, 0, 0, 1, 1, 2 };
                int[] Ji = { 0, 1, 2, 3, 5, 6, 0, 1, 2, 4, 0, 1, 2, 3, 0, 1, 2, 3, 1, 0, 1, 2, 3, 4, 12, 0, 1, 2, 0, 2, 2 };
                double[] ni = { 5.91599780322238E-05, -1.85465997137856E-03, 1.04190510480013E-02, 5.9864730203859E-03, -0.771391189901699, 1.72549765557036, -4.67076079846526E-04, 1.34533823384439E-02, -8.08094336805495E-02, 0.508139374365767, 1.28584643361683E-03, -1.63899353915435, 5.86938199318063, -2.92466667918613, -6.14076301499537E-03, 5.76199014049172, -12.1613320606788, 1.67637540957944, -7.44135838773463, 3.78168091437659E-02, 4.01432203027688, 16.0279837479185, 3.17848779347728, -3.58362310304853, -1159952.60446827, 0.199256573577909, -0.122270624794624, -19.1449143716586, -1.50448002905284E-02, 14.6407900162154, -3.2747778718823 };
                ps = p / 100;
                sigma = s / 5.3;
                omega = 0;
                for (i = 0; i <= 30; i++)
                {
                    omega = omega + ni[i] * Math.Pow(ps + 0.298, Ii[i]) * Math.Pow(sigma - 0.816, Ji[i]);
                }
                tempv3_ps = omega * 0.0088;
            }
            return tempv3_ps;
        }

        public double p3_hs(double h, double s)
        {
            double tempp3_hs = 0;
            //Supplementary Release on Backward Equations ( ) , p h s for Region 3,
            //Equations as a Function of h and s for the Region Boundaries, and an Equation
            //( ) sat , T hs for Region 4 of the IAPWS Industrial Formulation 1997 for the
            //Thermodynamic Properties of Water and Steam
            //2004
            //Section 3 Backward Functions p(h,s), T(h,s), and v(h,s) for Region 3
            int i = 0;

            double ps = 0;
            double sigma = 0;
            double eta = 0;
       //     double R = 0.461526;
       //     double tc = 647.096;
        //    double pc = 22.064;
        //    double rhoc = 322;
            if (s < 4.41202148223476)
            {
                //Subregion 3a
                //Eq 1, Table 3, Page 8
                int[] Ii = { 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 6, 7, 8, 10, 10, 14, 18, 20, 22, 22, 24, 28, 28, 32, 32 };
                int[] Ji = { 0, 1, 5, 0, 3, 4, 8, 14, 6, 16, 0, 2, 3, 0, 1, 4, 5, 28, 28, 24, 1, 32, 36, 22, 28, 36, 16, 28, 36, 16, 36, 10, 28 };
                double[] ni = { 7.70889828326934, -26.0835009128688, 267.416218930389, 17.2221089496844, -293.54233214597, 614.135601882478, -61056.2757725674, -65127225.1118219, 73591.9313521937, -11664650591.4191, 35.5267086434461, -596.144543825955, -475.842430145708, 69.6781965359503, 335.674250377312, 25052.6809130882, 146997.380630766, 5.38069315091534E+19, 1.43619827291346E+21, 3.64985866165994E+19, -2547.41561156775, 2.40120197096563E+27, -3.93847464679496E+29, 1.47073407024852E+24, -4.26391250432059E+31, 1.94509340621077E+38, 6.66212132114896E+23, 7.06777016552858E+33, 1.75563621975576E+41, 1.08408607429124E+28, 7.30872705175151E+43, 1.5914584739887E+24, 3.77121605943324E+40 };


                sigma = s / 4.4;
                eta = h / 2300;
                ps = 0;
                for (i = 0; i <= 32; i++)
                {
                    ps = ps + ni[i] * Math.Pow(eta - 1.01, Ii[i]) * Math.Pow(sigma - 0.75, Ji[i]);
                }
                tempp3_hs = ps * 99;
            }
            else
            {
                //Subregion 3b
                //Eq 2, Table 4, Page 8
                int[] Ii = { -12, -12, -12, -12, -12, -10, -10, -10, -10, -8, -8, -6, -6, -6, -6, -5, -4, -4, -4, -3, -3, -3, -3, -2, -2, -1, 0, 2, 2, 5, 6, 8, 10, 14, 14 };
                int[] Ji = { 2, 10, 12, 14, 20, 2, 10, 14, 18, 2, 8, 2, 6, 7, 8, 10, 4, 5, 8, 1, 3, 5, 6, 0, 1, 0, 3, 0, 1, 0, 1, 1, 1, 3, 7 };
                double[] ni = { 1.25244360717979E-13, -1.26599322553713E-02, 5.06878030140626, 31.7847171154202, -391041.161399932, -9.75733406392044E-11, -18.6312419488279, 510.973543414101, 373847.005822362, 2.99804024666572E-08, 20.0544393820342, -4.98030487662829E-06, -10.230180636003, 55.2819126990325, -206.211367510878, -7940.12232324823, 7.82248472028153, -58.6544326902468, 3550.73647696481, -1.15303107290162E-04, -1.75092403171802, 257.98168774816, -727.048374179467, 1.21644822609198E-04, 3.93137871762692E-02, 7.04181005909296E-03, -82.910820069811, -0.26517881813125, 13.7531682453991, -52.2394090753046, 2405.56298941048, -22736.1631268929, 89074.6343932567, -23923456.5822486, 5687958081.29714 };
                sigma = s / 5.3;
                eta = h / 2800;
                ps = 0;
                for (i = 0; i <= 34; i++)
                {
                    ps = ps + ni[i] * Math.Pow(eta - 0.681, Ii[i]) * Math.Pow(sigma - 0.792, Ji[i]);
                }
                tempp3_hs = 16.6 / ps;
            }
            return tempp3_hs;
        }

        public double h3_pT(double p, double T)
        {
            //Not avalible with IF 97
            //Solve function T3_ph-T=0 with half interval method.
            double Ts = 0;
            double Low_Bound = 0;
            double High_Bound = 0;
            double hs = 0;
            //ver2.6 Start corrected bug
            if (p < 22.06395) //Bellow tripple point
            {
                Ts = T4_p(p); //Saturation temperature
                if (T <= Ts) //Liquid side
                {
                    High_Bound = h4L_p(p); //Max h är liauid h.
                    Low_Bound = h1_pT(p, 623.15);
                }
                else
                {
                    Low_Bound = h4V_p(p); //Min h är Vapour h.
                    High_Bound = h2_pT(p, B23T_p(p));
                }
            }
            else //Above tripple point. R3 from R2 till R3.
            {
                Low_Bound = h1_pT(p, 623.15);
                High_Bound = h2_pT(p, B23T_p(p));
            }
            //ver2.6 End corrected bug
            Ts = T + 1;
            while (Math.Abs(T - Ts) > 0.000001)
            {
                hs = (Low_Bound + High_Bound) / 2;
                Ts = T3_ph(p, hs);
                if (Ts > T)
                {
                    High_Bound = hs;
                }
                else
                {
                    Low_Bound = hs;
                }
            }
            return hs;
        }

        public double T3_prho(double p, double rho)
        {
            //Solve by iteration. Observe that fo low temperatures this equation has 2 solutions.
            //Solve with half interval method
            double ps = 0;
            double Low_Bound = 0;
            double High_Bound = 0;
            double Ts = 0;
            Low_Bound = 623.15;
            High_Bound = 1073.15;
            while (Math.Abs(p - ps) > 0.00000001)
            {
                Ts = (Low_Bound + High_Bound) / 2;
                ps = p3_rhoT(rho, Ts);
                if (ps > p)
                {
                    High_Bound = Ts;
                }
                else
                {
                    Low_Bound = Ts;
                }
            }
            return Ts;
        }

        //***********************************************************************************************************
        //*2.4 Functions for region 4

        public double p4_T(double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //Section 8.1 The Saturation-Pressure Equation
            //Eq 30, Page 33
            double teta = 0;
            double a = 0;
            double b = 0;
            double c = 0;
            teta = T - 0.23855557567849 / (T - 650.17534844798);
            a = Math.Pow(teta, 2) + 1167.0521452767 * teta - 724213.16703206;
            b = -17.073846940092 * Math.Pow(teta, 2) + 12020.82470247 * teta - 3232555.0322333;
            c = 14.91510861353 * Math.Pow(teta, 2) - 4823.2657361591 * teta + 405113.40542057;
            return Math.Pow(2 * c / (-b + Math.Pow(Math.Pow(b, 2) - 4 * a * c, 0.5)), 4);
        }

        public double T4_p(double p)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //Section 8.2 The Saturation-Temperature Equation
            //Eq 31, Page 34
            double beta = 0;
            double e = 0;
            double f = 0;
            double g = 0;
            double d = 0;
            beta = Math.Pow(p, 0.25);
            e = Math.Pow(beta, 2) - 17.073846940092 * beta + 14.91510861353;
            f = 1167.0521452767 * Math.Pow(beta, 2) + 12020.82470247 * beta - 4823.2657361591;
            g = -724213.16703206 * Math.Pow(beta, 2) - 3232555.0322333 * beta + 405113.40542057;
            d = 2 * g / (-f - Math.Pow(Math.Pow(f, 2) - 4 * e * g, 0.5));
            return (650.17534844798 + d - Math.Pow(Math.Pow(650.17534844798 + d, 2) - 4 * (-0.23855557567849 + 650.17534844798 * d), 0.5)) / 2;
        }

        public double h4_s(double s)
        {
            double temph4_s = 0;
            //Supplementary Release on Backward Equations ( ) , p h s for Region 3,Equations as a Function of h and s for the Region Boundaries, and an Equation( ) sat , T hs for Region 4 of the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam
            //4 Equations for Region Boundaries Given Enthalpy and Entropy
            // Se picture page 14

            double eta = 0;
            double sigma = 0;
            double sigma1 = 0;
            double sigma2 = 0;
            int i = 0;
            if (s > -0.0001545495919 && s <= 3.77828134)
            {
                //hL1_s
                //Eq 3,Table 9,Page 16
                int[] Ii = { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 7, 8, 12, 12, 14, 14, 16, 20, 20, 22, 24, 28, 32, 32 };
                int[] Ji = { 14, 36, 3, 16, 0, 5, 4, 36, 4, 16, 24, 18, 24, 1, 4, 2, 4, 1, 22, 10, 12, 28, 8, 3, 0, 6, 8 };
                double[] ni = { 0.332171191705237, 6.11217706323496E-04, -8.82092478906822, -0.45562819254325, -2.63483840850452E-05, -22.3949661148062, -4.28398660164013, -0.616679338856916, -14.682303110404, 284.523138727299, -113.398503195444, 1156.71380760859, 395.551267359325, -1.54891257229285, 19.4486637751291, -3.57915139457043, -3.35369414148819, -0.66442679633246, 32332.1885383934, 3317.66744667084, -22350.1257931087, 5739538.75852936, 173.226193407919, -3.63968822121321E-02, 8.34596332878346E-07, 5.03611916682674, 65.5444787064505 };
                sigma = s / 3.8;
                eta = 0;
                for (i = 0; i <= 26; i++)
                {
                    eta = eta + ni[i] * Math.Pow(sigma - 1.09, Ii[i]) * Math.Pow(sigma + 0.0000366, Ji[i]);
                }
                temph4_s = eta * 1700;
            }
            else if (s > 3.77828134 && s <= 4.41202148223476)
            {
                //hL3_s
                //Eq 4,Table 10,Page 16
                int[] Ii = { 0, 0, 0, 0, 2, 3, 4, 4, 5, 5, 6, 7, 7, 7, 10, 10, 10, 32, 32 };
                int[] Ji = { 1, 4, 10, 16, 1, 36, 3, 16, 20, 36, 4, 2, 28, 32, 14, 32, 36, 0, 6 };
                double[] ni = { 0.822673364673336, 0.181977213534479, -0.011200026031362, -7.46778287048033E-04, -0.179046263257381, 4.24220110836657E-02, -0.341355823438768, -2.09881740853565, -8.22477343323596, -4.99684082076008, 0.191413958471069, 5.81062241093136E-02, -1655.05498701029, 1588.70443421201, -85.0623535172818, -31771.4386511207, -94589.0406632871, -1.3927384708869E-06, 0.63105253224098 };

                sigma = s / 3.8;
                eta = 0;
                for (i = 0; i <= 18; i++)
                {
                    eta = eta + ni[i] * Math.Pow(sigma - 1.09, Ii[i]) * Math.Pow(sigma + 0.0000366, Ji[i]);
                }
                temph4_s = eta * 1700;
            }
            else if (s > 4.41202148223476 && s <= 5.85)
            {
                //Section 4.4 Equations ( ) 2ab " h s and ( ) 2c3b "h s for the Saturated Vapor Line
                //Page 19, Eq 5
                //hV2c3b_s(s)
                int[] Ii = { 0, 0, 0, 1, 1, 5, 6, 7, 8, 8, 12, 16, 22, 22, 24, 36 };
                int[] Ji = { 0, 3, 4, 0, 12, 36, 12, 16, 2, 20, 32, 36, 2, 32, 7, 20 };
                double[] ni = { 1.04351280732769, -2.27807912708513, 1.80535256723202, 0.420440834792042, -105721.24483466, 4.36911607493884E+24, -328032702839.753, -6.7868676080427E+15, 7439.57464645363, -3.56896445355761E+19, 1.67590585186801E+31, -3.55028625419105E+37, 396611982166.538, -4.14716268484468E+40, 3.59080103867382E+18, -1.16994334851995E+40 };
                sigma = s / 5.9;
                eta = 0;
                for (i = 0; i <= 15; i++)
                {
                    eta = eta + ni[i] * Math.Pow(sigma - 1.02, Ii[i]) * Math.Pow(sigma - 0.726, Ji[i]);
                }
                temph4_s = Math.Pow(eta, 4) * 2800;
            }
            else if (s > 5.85 && s < 9.155759395)
            {
                //Section 4.4 Equations ( ) 2ab " h s and ( ) 2c3b "h s for the Saturated Vapor Line
                //Page 20, Eq 6
                int[] Ii = { 1, 1, 2, 2, 4, 4, 7, 8, 8, 10, 12, 12, 18, 20, 24, 28, 28, 28, 28, 28, 32, 32, 32, 32, 32, 36, 36, 36, 36, 36 };
                int[] Ji = { 8, 24, 4, 32, 1, 2, 7, 5, 12, 1, 0, 7, 10, 12, 32, 8, 12, 20, 22, 24, 2, 7, 12, 14, 24, 10, 12, 20, 22, 28 };
                double[] ni = { -524.581170928788, -9269472.18142218, -237.385107491666, 21077015581.2776, -23.9494562010986, 221.802480294197, -5104725.33393438, 1249813.96109147, 2000084369.96201, -815.158509791035, -157.612685637523, -11420042233.2791, 6.62364680776872E+15, -2.27622818296144E+18, -1.71048081348406E+31, 6.60788766938091E+15, 1.66320055886021E+22, -2.18003784381501E+29, -7.87276140295618E+29, 1.51062329700346E+31, 7957321.70300541, 1.31957647355347E+15, -3.2509706829914E+23, -4.18600611419248E+25, 2.97478906557467E+34, -9.53588761745473E+19, 1.66957699620939E+24, -1.75407764869978E+32, 3.47581490626396E+34, -7.10971318427851E+38 };
                sigma1 = s / 5.21;
                sigma2 = s / 9.2;
                eta = 0;
                for (i = 0; i <= 29; i++)
                {
                    eta = eta + ni[i] * Math.Pow(1 / sigma1 - 0.513, Ii[i]) * Math.Pow(sigma2 - 0.524, Ji[i]);
                }
                temph4_s = Math.Exp(eta) * 2800;
            }
            else
            {
                temph4_s = 9999999999;  // error // need to check how to handle this
            }
            return temph4_s;
        }

        public double p4_s(double s)
        {
            double tempp4_s = 0;
            //Uses h4_s and p_hs for the diffrent regions to determine p4_s
            double hsat;
            hsat = h4_s(s);
            if (s > -0.0001545495919 && s <= 3.77828134)
            {
                tempp4_s = p1_hs(hsat, s);
            }
            else if (s > 3.77828134 && s <= 5.210887663)
            {
                tempp4_s = p3_hs(hsat, s);
            }
            else if (s > 5.210887663 && s < 9.155759395)
            {
                tempp4_s = p2_hs(hsat, s);
            }
            else
            {
                tempp4_s = 9999999999;  // error // need to check how to handle this
            }
            return tempp4_s;
        }

        public double h4L_p(double p)
        {
            double temph4L_p = 0;
            double Low_Bound = 0;
            double High_Bound = 0;
            double hs = 0;
            double ps = 0;
            double Ts = 0;
            if (p > 0.000611657 && p < 22.06395)
            {
                Ts = T4_p(p);
                if (p < 16.529)
                {
                    temph4L_p = h1_pT(p, Ts);
                }
                else
                {
                    //Iterate to find the the backward solution of p3sat_h
                    Low_Bound = 1670.858218;
                    High_Bound = 2087.23500164864;
                    while (Math.Abs(p - ps) > 0.00001)
                    {
                        hs = (Low_Bound + High_Bound) / 2;
                        ps = p3sat_h(hs);
                        if (ps > p)
                        {
                            High_Bound = hs;
                        }
                        else
                        {
                            Low_Bound = hs;
                        }
                    }

                    temph4L_p = hs;
                }
            }
            else
            {
                temph4L_p = 9999999999;  // error // need to check how to handle this
            }
            return temph4L_p;
        }

        public double h4V_p(double p)
        {
            double temph4V_p = 0;
            double Low_Bound = 0;
            double High_Bound = 0;
            double hs = 0;
            double ps = 0;
            double Ts = 0;
            if (p > 0.000611657 && p < 22.06395)
            {
                Ts = T4_p(p);
                if (p < 16.529)
                {
                    temph4V_p = h2_pT(p, Ts);
                }
                else
                {
                    //Iterate to find the the backward solution of p3sat_h
                    Low_Bound = 2087.23500164864;
                    High_Bound = 2563.592004 + 5; //5 added to extrapolate to ensure even the border ==350°C solved.
                    while (Math.Abs(p - ps) > 0.000001)
                    {
                        hs = (Low_Bound + High_Bound) / 2;
                        ps = p3sat_h(hs);
                        if (ps < p)
                        {
                            High_Bound = hs;
                        }
                        else
                        {
                            Low_Bound = hs;
                        }
                    }
                    temph4V_p = hs;
                }
            }
            else
            {
                temph4V_p = 9999999999;  // error // need to check how to handle this
            }
            return temph4V_p;
        }

        public double x4_ph(double p, double h)
        {
            double tempx4_ph = 0;
            //Calculate vapour fraction from hL and hV for given p
            double h4v = 0;
            double h4l = 0;
            h4v = h4V_p(p);
            h4l = h4L_p(p);
            if (h > h4v)
            {
                tempx4_ph = 1;
            }
            else if (h < h4l)
            {
                tempx4_ph = 0;
            }
            else
            {
                tempx4_ph = (h - h4l) / (h4v - h4l);
            }
            return tempx4_ph;
        }

        public double x4_ps(double p, double s)
        {
            double tempx4_ps = 0;
            double ssV = 0;
            double ssL = 0;
            if (p < 16.529)
            {
                ssV = s2_pT(p, T4_p(p));
                ssL = s1_pT(p, T4_p(p));
            }
            else
            {
                ssV = s3_rhoT(1 / (v3_ph(p, h4V_p(p))), T4_p(p));
                ssL = s3_rhoT(1 / (v3_ph(p, h4L_p(p))), T4_p(p));
            }
            if (s < ssL)
            {
                tempx4_ps = 0;
            }
            else if (s > ssV)
            {
                tempx4_ps = 1;
            }
            else
            {
                tempx4_ps = (s - ssL) / (ssV - ssL);
            }
            return tempx4_ps;
        }

        public double T4_hs(double h, double s)
        {
            double tempT4_hs = 0;
            //Supplementary Release on Backward Equations ( ) , p h s for Region 3,
            //Chapter 5.3 page 30.
            //The if 97 function is only valid for part of region4. Use iteration outsida.

            double hL = 0;
            double Ts = 0;
            double ss = 0;
            double p = 0;
            double sigma = 0;
            double eta = 0;
            double teta = 0;
            double High_Bound = 0;
            double Low_Bound = 0;
            double PL = 0;
            double s4V = 0;
            double v4V = 0;
            double s4L = 0;
            double v4L = 0;
            double xs = 0;
            int i = 0;
            int[] Ii = { 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 8, 10, 10, 12, 14, 14, 16, 16, 18, 18, 18, 20, 28 };
            int[] Ji = { 0, 3, 12, 0, 1, 2, 5, 0, 5, 8, 0, 2, 3, 4, 0, 1, 1, 2, 4, 16, 6, 8, 22, 1, 20, 36, 24, 1, 28, 12, 32, 14, 22, 36, 24, 36 };
            double[] ni = { 0.179882673606601, -0.267507455199603, 1.162767226126, 0.147545428713616, -0.512871635973248, 0.421333567697984, 0.56374952218987, 0.429274443819153, -3.3570455214214, 10.8890916499278, -0.248483390456012, 0.30415322190639, -0.494819763939905, 1.07551674933261, 7.33888415457688E-02, 1.40170545411085E-02, -0.106110975998808, 1.68324361811875E-02, 1.25028363714877, 1013.16840309509, -1.51791558000712, 52.4277865990866, 23049.5545563912, 2.49459806365456E-02, 2107964.67412137, 366836848.613065, -144814105.365163, -1.7927637300359E-03, 4899556021.00459, 471.262212070518, -82929439019.8652, -1715.45662263191, 3557776.82973575, 586062760258.436, -12988763.5078195, 31724744937.1057 };


            if (s > 5.210887825 && s < 9.15546555571324)
            {
                sigma = s / 9.2;
                eta = h / 2800;
                teta = 0;
                for (i = 0; i <= 35; i++)
                {
                    teta = teta + ni[i] * Math.Pow(eta - 0.119, Ii[i]) * Math.Pow(sigma - 1.07, Ji[i]);
                }
                tempT4_hs = teta * 550;
            }
            else
            {
                //Function psat_h
                if (s > -0.0001545495919 && s <= 3.77828134)
                {
                    Low_Bound = 0.000611;
                    High_Bound = 165.291642526045;
                    while (Math.Abs(hL - h) > 0.00001 && Math.Abs(High_Bound - Low_Bound) > 0.0001)
                    {
                        PL = (Low_Bound + High_Bound) / 2;
                        Ts = T4_p(PL);
                        hL = h1_pT(PL, Ts);
                        if (hL > h)
                        {
                            High_Bound = PL;
                        }
                        else
                        {
                            Low_Bound = PL;
                        }
                    }
                }
                else if (s > 3.77828134 && s <= 4.41202148223476)
                {
                    PL = p3sat_h(h);
                }
                else if (s > 4.41202148223476 && s <= 5.210887663)
                {
                    PL = p3sat_h(h);
                }
                Low_Bound = 0.000611;
                High_Bound = PL;
                while (Math.Abs(s - ss) > 0.000001 && Math.Abs(High_Bound - Low_Bound) > 0.0000001)
                {
                    p = (Low_Bound + High_Bound) / 2;

                    //Calculate s4_ph
                    Ts = T4_p(p);
                    xs = x4_ph(p, h);
                    if (p < 16.529)
                    {
                        s4V = s2_pT(p, Ts);
                        s4L = s1_pT(p, Ts);
                    }
                    else
                    {
                        v4V = v3_ph(p, h4V_p(p));
                        s4V = s3_rhoT(1 / v4V, Ts);
                        v4L = v3_ph(p, h4L_p(p));
                        s4L = s3_rhoT(1 / v4L, Ts);
                    }
                    ss = (xs * s4V + (1 - xs) * s4L);

                    if (ss < s)
                    {
                        High_Bound = p;
                    }
                    else
                    {
                        Low_Bound = p;
                    }
                }
                tempT4_hs = T4_p(p);
            }
            return tempT4_hs;
        }


        //***********************************************************************************************************
        //*2.5 Functions for region 5

        public double h5_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //Basic Equation for Region 5
            //Eq 32,33, Page 36, Tables 37-41





            double tau = 0;
            double gamma0_tau = 0;
            double gammar_tau = 0;
            int i = 0;
            double R = 0.461526; //kJ/(kg K)
            int[] Ji0 = { 0, 1, -3, -2, -1, 2 };
            double[] ni0 = { -13.179983674201, 6.8540841634434, -0.024805148933466, 0.36901534980333, -3.1161318213925, -0.32961626538917 };
            int[] Iir = { 1, 1, 1, 2, 3 };
            int[] Jir = { 0, 1, 3, 9, 3 };
            double[] nir = { -1.2563183589592E-04, 2.1774678714571E-03, -0.004594282089991, -3.9724828359569E-06, 1.2919228289784E-07 };
            tau = 1000 / T;
            gamma0_tau = 0;
            for (i = 0; i <= 5; i++)
            {
                gamma0_tau = gamma0_tau + ni0[i] * Ji0[i] * Math.Pow(tau, Ji0[i] - 1);
            }
            gammar_tau = 0;
            for (i = 0; i <= 4; i++)
            {
                gammar_tau = gammar_tau + nir[i] * Math.Pow(p, Iir[i]) * Jir[i] * Math.Pow(tau, Jir[i] - 1);
            }
            return R * T * tau * (gamma0_tau + gammar_tau);
        }

        public double v5_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //Basic Equation for Region 5
            //Eq 32,33, Page 36, Tables 37-41

            double tau = 0;
            double gamma0_pi = 0;
            double gammar_pi = 0;
            int i = 0;
            double R = 0.461526; //kJ/(kg K)
            int[] Ji0 = { 0, 1, -3, -2, -1, 2 };
            double[] ni0 = { -13.179983674201, 6.8540841634434, -0.024805148933466, 0.36901534980333, -3.1161318213925, -0.32961626538917 };
            int[] Iir = { 1, 1, 1, 2, 3 };
            int[] Jir = { 0, 1, 3, 9, 3 };
            double[] nir = { -1.2563183589592E-04, 2.1774678714571E-03, -0.004594282089991, -3.9724828359569E-06, 1.2919228289784E-07 };
            tau = 1000 / T;
            gamma0_pi = 1 / p;
            gammar_pi = 0;
            for (i = 0; i <= 4; i++)
            {
                gammar_pi = gammar_pi + nir[i] * Iir[i] * Math.Pow(p, Iir[i] - 1) * Math.Pow(tau, Jir[i]);
            }
            return R * T / p * p * (gamma0_pi + gammar_pi) / 1000;
        }

        public double u5_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //Basic Equation for Region 5
            //Eq 32,33, Page 36, Tables 37-41

            double tau = 0;
            double gamma0_pi = 0;
            double gammar_pi = 0;
            double gamma0_tau = 0;
            double gammar_tau = 0;
            int i = 0;
            double R = 0.461526; //kJ/(kg K)
            int[] Ji0 = { 0, 1, -3, -2, -1, 2 };
            double[] ni0 = { -13.179983674201, 6.8540841634434, -0.024805148933466, 0.36901534980333, -3.1161318213925, -0.32961626538917 };
            int[] Iir = { 1, 1, 1, 2, 3 };
            int[] Jir = { 0, 1, 3, 9, 3 };
            double[] nir = { -1.2563183589592E-04, 2.1774678714571E-03, -0.004594282089991, -3.9724828359569E-06, 1.2919228289784E-07 };
            tau = 1000 / T;
            gamma0_pi = 1 / p;
            gamma0_tau = 0;
            for (i = 0; i <= 5; i++)
            {
                gamma0_tau = gamma0_tau + ni0[i] * Ji0[i] * Math.Pow(tau, Ji0[i] - 1);
            }
            gammar_pi = 0;
            gammar_tau = 0;
            for (i = 0; i <= 4; i++)
            {
                gammar_pi = gammar_pi + nir[i] * Iir[i] * Math.Pow(p, Iir[i] - 1) * Math.Pow(tau, Jir[i]);
                gammar_tau = gammar_tau + nir[i] * Math.Pow(p, Iir[i]) * Jir[i] * Math.Pow(tau, Jir[i] - 1);
            }
            return R * T * (tau * (gamma0_tau + gammar_tau) - p * (gamma0_pi + gammar_pi));
        }

        public double Cp5_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //Basic Equation for Region 5
            //Eq 32,33, Page 36, Tables 37-41

            double tau = 0;
            double gamma0_tautau = 0;
            double gammar_tautau = 0;
            int i = 0;
            double R = 0.461526; //kJ/(kg K)
            int[] Ji0 = { 0, 1, -3, -2, -1, 2 };
            double[] ni0 = { -13.179983674201, 6.8540841634434, -0.024805148933466, 0.36901534980333, -3.1161318213925, -0.32961626538917 };
            int[] Iir = { 1, 1, 1, 2, 3 };
            int[] Jir = { 0, 1, 3, 9, 3 };
            double[] nir = { -1.2563183589592E-04, 2.1774678714571E-03, -0.004594282089991, -3.9724828359569E-06, 1.2919228289784E-07 };
            tau = 1000 / T;
            gamma0_tautau = 0;
            for (i = 0; i <= 5; i++)
            {
                gamma0_tautau = gamma0_tautau + ni0[i] * Ji0[i] * (Ji0[i] - 1) * Math.Pow(tau, Ji0[i] - 2);
            }
            gammar_tautau = 0;
            for (i = 0; i <= 4; i++)
            {
                gammar_tautau = gammar_tautau + nir[i] * Math.Pow(p, Iir[i]) * Jir[i] * (Jir[i] - 1) * Math.Pow(tau, Jir[i] - 2);
            }
            return -R * Math.Pow(tau, 2) * (gamma0_tautau + gammar_tautau);
        }

        public double s5_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //Basic Equation for Region 5
            //Eq 32,33, Page 36, Tables 37-41

            double tau = 0;
            double gamma0 = 0;
            double gamma0_tau = 0;
            double gammar = 0;
            double gammar_tau = 0;
            int i = 0;
            double R = 0.461526; //kJ/(kg K)
            int[] Ji0 = { 0, 1, -3, -2, -1, 2 };
            double[] ni0 = { -13.179983674201, 6.8540841634434, -0.024805148933466, 0.36901534980333, -3.1161318213925, -0.32961626538917 };
            int[] Iir = { 1, 1, 1, 2, 3 };
            int[] Jir = { 0, 1, 3, 9, 3 };
            double[] nir = { -1.2563183589592E-04, 2.1774678714571E-03, -0.004594282089991, -3.9724828359569E-06, 1.2919228289784E-07 };

            tau = 1000 / T;
            gamma0 = Math.Log(p);
            gamma0_tau = 0;
            for (i = 0; i <= 5; i++)
            {
                gamma0_tau = gamma0_tau + ni0[i] * Ji0[i] * Math.Pow(tau, Ji0[i] - 1);
                gamma0 = gamma0 + ni0[i] * Math.Pow(tau, Ji0[i]);
            }
            gammar = 0;
            gammar_tau = 0;
            for (i = 0; i <= 4; i++)
            {
                gammar = gammar + nir[i] * Math.Pow(p, Iir[i]) * Math.Pow(tau, Jir[i]);
                gammar_tau = gammar_tau + nir[i] * Math.Pow(p, Iir[i]) * Jir[i] * Math.Pow(tau, Jir[i] - 1);
            }
            return R * (tau * (gamma0_tau + gammar_tau) - (gamma0 + gammar));
        }

        public double Cv5_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //Basic Equation for Region 5
            //Eq 32,33, Page 36, Tables 37-41

            double tau = 0;
            double gamma0_tautau = 0;
            double gammar_pi = 0;
            double gammar_pitau = 0;
            double gammar_pipi = 0;
            double gammar_tautau = 0;
            int i = 0;
            double R = 0.461526; //kJ/(kg K)
            int[] Ji0 = { 0, 1, -3, -2, -1, 2 };
            double[] ni0 = { -13.179983674201, 6.8540841634434, -0.024805148933466, 0.36901534980333, -3.1161318213925, -0.32961626538917 };
            int[] Iir = { 1, 1, 1, 2, 3 };
            int[] Jir = { 0, 1, 3, 9, 3 };
            double[] nir = { -1.2563183589592E-04, 2.1774678714571E-03, -0.004594282089991, -3.9724828359569E-06, 1.2919228289784E-07 };

            tau = 1000 / T;
            gamma0_tautau = 0;
            for (i = 0; i <= 5; i++)
            {
                gamma0_tautau = gamma0_tautau + ni0[i] * (Ji0[i] - 1) * Ji0[i] * Math.Pow(tau, Ji0[i] - 2);
            }
            gammar_pi = 0;
            gammar_pitau = 0;
            gammar_pipi = 0;
            gammar_tautau = 0;
            for (i = 0; i <= 4; i++)
            {
                gammar_pi = gammar_pi + nir[i] * Iir[i] * Math.Pow(p, Iir[i] - 1) * Math.Pow(tau, Jir[i]);
                gammar_pitau = gammar_pitau + nir[i] * Iir[i] * Math.Pow(p, Iir[i] - 1) * Jir[i] * Math.Pow(tau, Jir[i] - 1);
                gammar_pipi = gammar_pipi + nir[i] * Iir[i] * (Iir[i] - 1) * Math.Pow(p, Iir[i] - 2) * Math.Pow(tau, Jir[i]);
                gammar_tautau = gammar_tautau + nir[i] * Math.Pow(p, Iir[i]) * Jir[i] * (Jir[i] - 1) * Math.Pow(tau, Jir[i] - 2);
            }
            return R * (-(Math.Pow(tau, 2) * (gamma0_tautau + gammar_tautau)) - Math.Pow(1 + p * gammar_pi - tau * p * gammar_pitau, 2) / (1 - Math.Pow(p, 2) * gammar_pipi));

        }

        public double w5_pT(double p, double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
            //Basic Equation for Region 5
            //Eq 32,33, Page 36, Tables 37-41

            double tau = 0;
            double gamma0_tautau = 0;
            double gammar_pi = 0;
            double gammar_pitau = 0;
            double gammar_pipi = 0;
            double gammar_tautau = 0;
            int i = 0;
            double R = 0.461526; //kJ/(kg K)
            int[] Ji0 = { 0, 1, -3, -2, -1, 2 };
            double[] ni0 = { -13.179983674201, 6.8540841634434, -0.024805148933466, 0.36901534980333, -3.1161318213925, -0.32961626538917 };
            int[] Iir = { 1, 1, 1, 2, 3 };
            int[] Jir = { 0, 1, 3, 9, 3 };
            double[] nir = { -1.2563183589592E-04, 2.1774678714571E-03, -0.004594282089991, -3.9724828359569E-06, 1.2919228289784E-07 };

            tau = 1000 / T;
            gamma0_tautau = 0;
            for (i = 0; i <= 5; i++)
            {
                gamma0_tautau = gamma0_tautau + ni0[i] * (Ji0[i] - 1) * Ji0[i] * Math.Pow(tau, Ji0[i] - 2);
            }
            gammar_pi = 0;
            gammar_pitau = 0;
            gammar_pipi = 0;
            gammar_tautau = 0;
            for (i = 0; i <= 4; i++)
            {
                gammar_pi = gammar_pi + nir[i] * Iir[i] * Math.Pow(p, Iir[i] - 1) * Math.Pow(tau, Jir[i]);
                gammar_pitau = gammar_pitau + nir[i] * Iir[i] * Math.Pow(p, Iir[i] - 1) * Jir[i] * Math.Pow(tau, Jir[i] - 1);
                gammar_pipi = gammar_pipi + nir[i] * Iir[i] * (Iir[i] - 1) * Math.Pow(p, Iir[i] - 2) * Math.Pow(tau, Jir[i]);
                gammar_tautau = gammar_tautau + nir[i] * Math.Pow(p, Iir[i]) * Jir[i] * (Jir[i] - 1) * Math.Pow(tau, Jir[i] - 2);
            }
            return Math.Pow(1000 * R * T * (1 + 2 * p * gammar_pi + Math.Pow(p, 2) * Math.Pow(gammar_pi, 2)) / ((1 - Math.Pow(p, 2) * gammar_pipi) + Math.Pow(1 + p * gammar_pi - tau * p * gammar_pitau, 2) / (Math.Pow(tau, 2) * (gamma0_tautau + gammar_tautau))), 0.5);
        }


        public double T5_ph(double p, double h)
        {
            //Solve with half interval method
            double Low_Bound = 0;
            double High_Bound = 0;
            double Ts = 0;
            double hs = 0;
            Low_Bound = 1073.15;
            High_Bound = 2273.15;
            while (Math.Abs(h - hs) > 0.00001)
            {
                Ts = (Low_Bound + High_Bound) / 2;
                hs = h5_pT(p, Ts);
                if (hs > h)
                {
                    High_Bound = Ts;
                }
                else
                {
                    Low_Bound = Ts;
                }
            }
            return Ts;
        }

        public double T5_ps(double p, double s)
        {
            //Solve with half interval method
            double Low_Bound = 0;
            double High_Bound = 0;
            double Ts = 0;
            double ss = 0;
            Low_Bound = 1073.15;
            High_Bound = 2273.15;
            while (Math.Abs(s - ss) > 0.00001)
            {
                Ts = (Low_Bound + High_Bound) / 2;
                ss = s5_pT(p, Ts);
                if (ss > s)
                {
                    High_Bound = Ts;
                }
                else
                {
                    Low_Bound = Ts;
                }
            }
            return Ts;
        }

        public double T5_prho(double p, double rho)
        {
            //Solve by iteration. Observe that fo low temperatures this equation has 2 solutions.
            //Solve with half interval method
            double Low_Bound = 0;
            double High_Bound = 0;
            double Ts = 0;
            double rhos = 0;
            Low_Bound = 1073.15;
            High_Bound = 2073.15;
            while (Math.Abs(rho - rhos) > 0.000001)
            {
                Ts = (Low_Bound + High_Bound) / 2;
                rhos = 1 / v2_pT(p, Ts);
                if (rhos < rho)
                {
                    High_Bound = Ts;
                }
                else
                {
                    Low_Bound = Ts;
                }
            }
            return Ts;
        }

        //***********************************************************************************************************
        //*3 Region Selection
        //***********************************************************************************************************
        //*3.1 Regions as a function of pT


        public int region_pT(double p, double T)
        {
            int tempregion_pT = 0;
            double ps = 0;
            if (T > 1073.15 && p < 10 && T < 2273.15 && p > 0.000611)
            {
                tempregion_pT = 5;
            }
            else if (T <= 1073.15 && T > 273.15 && p <= 100 && p > 0.000611)
            {
                if (T > 623.15)
                {
                    if (p > B23p_T(T))
                    {
                        tempregion_pT = 3;
                        if (T < 647.096)
                        {
                            ps = p4_T(T);
                            if (Math.Abs(p - ps) < 0.00001)
                            {
                                tempregion_pT = 4;
                            }
                        }
                    }
                    else
                    {
                        tempregion_pT = 2;
                    }
                }
                else
                {
                    ps = p4_T(T);
                    if (Math.Abs(p - ps) < 0.00001)
                    {
                        tempregion_pT = 4;
                    }
                    else if (p > ps)
                    {
                        tempregion_pT = 1;
                    }
                    else
                    {
                        tempregion_pT = 2;
                    }
                }
            }
            else
            {
                tempregion_pT = 0; //**Error, Outside valid area
            }
            return tempregion_pT;
        }


        //***********************************************************************************************************
        //*3.2 Regions as a function of ph

        public int region_ph(double p, double h)
        {
            double hL = 0;
            double hV = 0;
            double h_45 = 0;
            double h_5u = 0;
            double Ts = 0;
            //Check if outside pressure limits
            if ((double)(p) < 0.000611657 || (int)(p) > 100)
            {
                return 0;
            }

            //Check if outside low h.
            if ((double)(h) < 0.963 * p + 2.2) //Linear adaption to h1_pt()+2 to speed up calcualations.
            {
                if (h < h1_pT(p, 273.15))
                {
                    return 0;
                }
            }

            if ((double)(p) < 16.5292) //Bellow region 3,Check region 1,4,2,5
            {
                //Check Region 1
                Ts = T4_p(p);
                hL = 109.6635 * Math.Log(p) + 40.3481 * p + 734.58; //Approximate function for hL_p
                if (Math.Abs(h - hL) < 100) //If approximate is not god enough use real function
                {
                    hL = h1_pT(p, Ts);
                }
                if ((double)(h) <= hL)
                {
                    return 1;
                }
                //Check Region 4
                hV = 45.1768 * Math.Log(p) - 20.158 * p + 2804.4; //Approximate function for hV_p
                if (Math.Abs(h - hV) < 50) //If approximate is not god enough use real function
                {
                    hV = h2_pT(p, Ts);
                }
                if ((double)(h) < hV)
                {
                    return 4;
                }
                //Check upper limit of region 2 Quick Test
                if ((int)(h) < 4000)
                {
                    return 2;
                }
                //Check region 2 (Real value)
                h_45 = h2_pT(p, 1073.15);
                if ((double)(h) <= h_45)
                {
                    return 2;
                }
                //Check region 5
                if ((int)(p) > 10)
                {
                    return 0;
                }
                h_5u = h5_pT(p, 2273.15);
                if ((double)(h) < h_5u)
                {
                    return 5;
                }
                return 0;
            }
            else //For p>16.5292
            {
                //Check if in region1
                if (h < h1_pT(p, 623.15))
                {
                    return 1;
                }
                //Check if in region 3 or 4 (Bellow Reg 2)
                if (h < h2_pT(p, B23T_p(p)))
                {
                    //Region 3 or 4
                    if (p > p3sat_h(h))
                    {
                        return 3;
                    }
                    else
                    {
                        return 4;
                    }
                }
                //Check if region 2
                if (h < h2_pT(p, 1073.15))
                {
                    return 2;
                }
            }
            return 0;
        }


        //***********************************************************************************************************
        //*3.3 Regions as a function of ps

        public int region_ps(double p, double s)
        {
            double ss = 0;
            if (p < 0.000611657 || p > 100 || (int)(s) < 0 || s > s5_pT(p, 2273.15))
            {
                return 0;
            }

            //Check region 5
            if (s > s2_pT(p, 1073.15))
            {
                if (p <= 10)
                {
                    return 5;
                }
                else
                {
                    return 0;
                }
            }

            //Check region 2
            if (p > 16.529)
            {
                ss = s2_pT(p, B23T_p(p)); //Between 5.047 and 5.261. Use to speed up!
            }
            else
            {
                ss = s2_pT(p, T4_p(p));
            }
            if ((double)(s) > ss)
            {
                return 2;
            }

            //Check region 3
            ss = s1_pT(p, 623.15);
            if (p > 16.529 && (double)(s) > ss)
            {
                if (p > p3sat_s(s))
                {
                    return 3;
                }
                else
                {
                    return 4;
                }
            }

            //Check region 4 (Not inside region 3)
            if (p < 16.529 && s > s1_pT(p, T4_p(p)))
            {
                return 4;
            }

            //Check region 1
            if (p > 0.000611657 && s > s1_pT(p, 273.15))
            {
                return 1;
            }
            return 1;
        }


        //***********************************************************************************************************
        //*3.4 Regions as a function of hs

        public int Region_hs(double h, double s)
        {
            double TMax = 0;
            double hMax = 0;
            double hB = 0;
            double hL = 0;
            double hV = 0;
            double vmax = 0;
            double Tmin = 0;
            double hMin = 0;
            if ((double)(s) < -0.0001545495919)
            {
                return 0;
            }
            //Check linear adaption to p=0.000611. If bellow region 4.
            hMin = (((-0.0415878 - 2500.89262) / (-0.00015455 - 9.155759)) * s);
            if ((double)(s) < 9.155759395 && (double)(h) < hMin)
            {
                return 0;
            }

            //******Kolla 1 eller 4. (+liten bit över B13)
            if ((double)(s) >= -0.0001545495919 && (double)(s) <= 3.77828134)
            {
                if (h < h4_s(s))
                {
                    return 4;
                }
                else if ((double)(s) < 3.397782955) //100MPa line is limiting
                {
                    TMax = T1_ps(100, s);
                    hMax = h1_pT(100, TMax);
                    if ((double)(h) < hMax)
                    {
                        return 1;
                    }
                    else
                    {
                        return 0;
                    }
                }
                else //The point is either in region 4,1,3. Check B23
                {
                    hB = hB13_s(s);
                    if ((double)(h) < hB)
                    {
                        return 1;
                    }
                    TMax = T3_ps(100, s);
                    vmax = v3_ps(100, s);
                    hMax = h3_rhoT(1 / vmax, TMax);
                    if ((double)(h) < hMax)
                    {
                        return 3;
                    }
                    else
                    {
                        return 0;
                    }
                }
            }

            //******Kolla region 2 eller 4. (Övre delen av område b23-> max)
            if ((double)(s) >= 5.260578707 && (double)(s) <= 11.9212156897728)
            {
                if ((double)(s) > 9.155759395) //Above region 4
                {
                    Tmin = T2_ps(0.000611, s);
                    hMin = h2_pT(0.000611, Tmin);
                    //Function adapted to h(1073.15,s)
                    hMax = -0.07554022 * Math.Pow(s, 4) + 3.341571 * Math.Pow(s, 3) - 55.42151 * Math.Pow(s, 2) + 408.515 * s + 3031.338;
                    if ((double)(h) > hMin && (double)(h) < hMax)
                    {
                        return 2;
                    }
                    else
                    {
                        return 0;
                    }
                }


                hV = h4_s(s);

                if ((double)(h) < hV) //Region 4. Under region 3.
                {
                    return 4;
                }
                if ((double)(s) < 6.04048367171238)
                {
                    TMax = T2_ps(100, s);
                    hMax = h2_pT(100, TMax);
                }
                else
                {
                    //Function adapted to h(1073.15,s)
                    hMax = -2.988734 * Math.Pow(s, 4) + 121.4015 * Math.Pow(s, 3) - 1805.15 * Math.Pow(s, 2) + 11720.16 * s - 23998.33;
                }
                if ((double)(h) < hMax) //Region 2. Över region 4.
                {
                    return 2;
                }
                else
                {
                    return 0;
                }
            }

            //Kolla region 3 eller 4. Under kritiska punkten.
            if ((double)(s) >= 3.77828134 && (double)(s) <= 4.41202148223476)
            {
                hL = h4_s(s);
                if ((double)(h) < hL)
                {
                    return 4;
                }
                TMax = T3_ps(100, s);
                vmax = v3_ps(100, s);
                hMax = h3_rhoT(1 / vmax, TMax);
                if ((double)(h) < hMax)
                {
                    return 3;
                }
                else
                {
                    return 0;
                }
            }

            //Kolla region 3 eller 4 från kritiska punkten till övre delen av b23
            if ((double)(s) >= 4.41202148223476 && (double)(s) <= 5.260578707)
            {
                hV = h4_s(s);
                if ((double)(h) < hV)
                {
                    return 4;
                }
                //Kolla om vi är under b23 giltighetsområde.
                if ((double)(s) <= 5.048096828)
                {
                    TMax = T3_ps(100, s);
                    vmax = v3_ps(100, s);
                    hMax = h3_rhoT(1 / vmax, TMax);
                    if ((double)(h) < hMax)
                    {
                        return 3;
                    }
                    else
                    {
                        return 0;
                    }
                }
                else //Inom området för B23 i s led.
                {
                    if ((double)(h) > 2812.942061) //Ovanför B23 i h_led
                    {
                        if ((double)(s) > 5.09796573397125)
                        {
                            TMax = T2_ps(100, s);
                            hMax = h2_pT(100, TMax);
                            if ((double)(h) < hMax)
                            {
                                return 2;
                            }
                            else
                            {
                                return 0;
                            }
                        }
                        else
                        {
                            return 0;
                        }
                    }
                    if ((double)(h) < 2563.592004) //Nedanför B23 i h_led men vi har redan kollat ovanför hV2c3b
                    {
                        return 3;
                    }
                    //Vi är inom b23 området i både s och h led.
                    if (p2_hs(h, s) > B23p_T(TB23_hs(h, s)))
                    {
                        return 3;
                    }
                    else
                    {
                        return 2;
                    }
                }
            }
            return 9;  // error // need to check how to handle this
        }


        //***********************************************************************************************************
        //*3.5 Regions as a function of p and rho

        public int Region_prho(double p, double rho)
        {
            double v;
            v = 1 / rho;
            if (p < 0.000611657 || p > 100)
            {
                return 0;
            }
            if (p < 16.5292) //Bellow region 3, Check region 1,4,2
            {
                if (v < v1_pT(p, 273.15)) //Observe that this is not actually min of v. Not valid Water of 4°C is ligther.
                {
                    return 0;
                }
                if (v <= v1_pT(p, T4_p(p)))
                {
                    return 1;
                }
                if (v < v2_pT(p, T4_p(p)))
                {
                    return 4;
                }
                if (v <= v2_pT(p, 1073.15))
                {
                    return 2;
                }
                if (p > 10) //Above region 5
                {
                    return 0;
                }
                if (v <= v5_pT(p, 2073.15))
                {
                    return 5;
                }
            }
            else //Check region 1,3,4,3,2 (Above the lowest point of region 3.)
            {
                if (v < v1_pT(p, 273.15)) //Observe that this is not actually min of v. Not valid Water of 4°C is ligther.
                {
                    return 0;
                }
                if (v < v1_pT(p, 623.15))
                {
                    return 1;
                }
                //Check if in region 3 or 4 (Bellow Reg 2)
                if (v < v2_pT(p, B23T_p(p)))
                {
                    //Region 3 or 4
                    if (p > 22.064) //Above region 4
                    {
                        return 3;
                    }
                    if (v < v3_ph(p, h4L_p(p)) || v > v3_ph(p, h4V_p(p))) //Uses iteration!!
                    {
                        return 3;
                    }
                    else
                    {
                        return 4;
                    }
                }
                //Check if region 2
                if (v < v2_pT(p, 1073.15))
                {
                    return 2;
                }
            }

            return 0;
        }


        //***********************************************************************************************************
        //*4 Region Borders
        //***********************************************************************************************************
        //***********************************************************************************************************
        //*4.1 Boundary between region 2 and 3.


        public double B23p_T(double T)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam
            //1997
            //Section 4 Auxiliary Equation for the Boundary between Regions 2 and 3
            //Eq 5, Page 5
            return 348.05185628969 - 1.1671859879975 * T + 1.0192970039326E-03 * Math.Pow(T, 2);
        }

        public double B23T_p(double p)
        {
            //Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam
            //1997
            //Section 4 Auxiliary Equation for the Boundary between Regions 2 and 3
            //Eq 6, Page 6
            return 572.54459862746 + Math.Pow((p - 13.91883977887) / 1.0192970039326E-03, 0.5);
        }


        //***********************************************************************************************************
        //*4.2 Region 3. pSat_h and pSat_s

        public double p3sat_h(double h)
        {
            //Revised Supplementary Release on Backward Equations for the Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam
            //2004
            //Section 4 Boundary Equations psat(h) and psat(s) for the Saturation Lines of Region 3
            //Se pictures Page 17, Eq 10, Table 17, Page 18

            double ps = 0;
            int i = 0;
            int[] Ii = { 0, 1, 1, 1, 1, 5, 7, 8, 14, 20, 22, 24, 28, 36 };
            int[] Ji = { 0, 1, 3, 4, 36, 3, 0, 24, 16, 16, 3, 18, 8, 24 };
            double[] ni = { 0.600073641753024, -9.36203654849857, 24.6590798594147, -107.014222858224, -91582131580576.8, -8623.32011700662, -23.5837344740032, 2.52304969384128E+17, -3.89718771997719E+18, -3.33775713645296E+22, 35649946963.6328, -1.48547544720641E+26, 3.30611514838798E+18, 8.13641294467829E+37 };
            h = h / 2600;
            ps = 0;
            for (i = 0; i <= 13; i++)
            {
                ps = ps + ni[i] * Math.Pow(h - 1.02, Ii[i]) * Math.Pow(h - 0.608, Ji[i]);
            }
            return ps * 22;
        }

        public double p3sat_s(double s)
        {

            double sigma = 0;
            double p = 0;
            int i = 0;
            int[] Ii = { 0, 1, 1, 4, 12, 12, 16, 24, 28, 32 };
            int[] Ji = { 0, 1, 32, 7, 4, 14, 36, 10, 0, 18 };
            double[] ni = { 0.639767553612785, -12.9727445396014, -2.24595125848403E+15, 1774667.41801846, 7170793495.71538, -3.78829107169011E+17, -9.55586736431328E+34, 1.87269814676188E+23, 119254746466.473, 1.10649277244882E+36 };

            sigma = s / 5.2;
            p = 0;
            for (i = 0; i <= 9; i++)
            {
                p = p + ni[i] * Math.Pow(sigma - 1.03, Ii[i]) * Math.Pow(sigma - 0.699, Ji[i]);
            }
            return p * 22;
        }

        //***********************************************************************************************************
        //4.3 Region boundary 1to3 and 3to2 as a functions of s

        public double hB13_s(double s)
        {
            //Supplementary Release on Backward Equations ( ) , p h s for Region 3,
            //Chapter 4.5 page 23.

            double sigma = 0;
            double eta = 0;
            int i = 0;
            int[] Ii = { 0, 1, 1, 3, 5, 6 };
            int[] Ji = { 0, -2, 2, -12, -4, -3 };
            double[] ni = { 0.913965547600543, -4.30944856041991E-05, 60.3235694765419, 1.17518273082168E-18, 0.220000904781292, -69.0815545851641 };

            eta = 0;
            for (i = 0; i <= 5; i++)
            {
                eta = eta + ni[i] * Math.Pow(sigma - 0.884, Ii[i]) * Math.Pow(sigma - 0.864, Ji[i]);
            }
            return eta * 1700;
        }

        public double TB23_hs(double h, double s)
        {
            //Supplementary Release on Backward Equations ( ) , p h s for Region 3,
            //Chapter 4.6 page 25.

            double sigma = 0;
            double eta = 0;
            double teta = 0;
            int i = 0;
            int[] Ii = { -12, -10, -8, -4, -3, -2, -2, -2, -2, 0, 1, 1, 1, 3, 3, 5, 6, 6, 8, 8, 8, 12, 12, 14, 14 };
            int[] Ji = { 10, 8, 3, 4, 3, -6, 2, 3, 4, 0, -3, -2, 10, -2, -1, -5, -6, -3, -8, -2, -1, -12, -1, -12, 1 };
            double[] ni = { 6.2909626082981E-04, -8.23453502583165E-04, 5.15446951519474E-08, -1.17565945784945, 3.48519684726192, -5.07837382408313E-12, -2.84637670005479, -2.36092263939673, 6.01492324973779, 1.48039650824546, 3.60075182221907E-04, -1.26700045009952E-02, -1221843.32521413, 0.149276502463272, 0.698733471798484, -2.52207040114321E-02, 1.47151930985213E-02, -1.08618917681849, -9.36875039816322E-04, 81.9877897570217, -182.041861521835, 2.61907376402688E-06, -29162.6417025961, 1.40660774926165E-05, 7832370.62349385 };
            sigma = s / 5.3;
            eta = h / 3000;
            teta = 0;
            for (i = 0; i <= 24; i++)
            {
                teta = teta + ni[i] * Math.Pow(eta - 0.727, Ii[i]) * Math.Pow(sigma - 0.864, Ji[i]);
            }
            return teta * 900;
        }


        //***********************************************************************************************************
        //*5 Transport properties
        //***********************************************************************************************************
        //*5.1 Viscosity (IAPWS formulation 1985, Revised 2003)
        //***********************************************************************************************************


        public double my_AllRegions_pT(double p, double T)
        {
            double rho = 0;
            double Ts = 0;
            double ps = 0;
            double my0 = 0;
            double sum = 0;
            double my1 = 0;
            double rhos = 0;
            int i = 0;
            double[] h0 = { 0.5132047, 0.3205656, 0, 0, -0.7782567, 0.1885447 };
            double[] h1 = { 0.2151778, 0.7317883, 1.241044, 1.476783, 0, 0 };
            double[] h2 = { -0.2818107, -1.070786, -1.263184, 0, 0, 0 };
            double[] h3 = { 0.1778064, 0.460504, 0.2340379, -0.4924179, 0, 0 };
            double[] h4 = { -0.0417661, 0, 0, 0.1600435, 0, 0 };
            double[] h5 = { 0, -0.01578386, 0, 0, 0, 0 };
            double[] h6 = { 0, 0, 0, -0.003629481, 0, 0 };


            //Calcualte density.
            switch (region_pT(p, T))
            {
                case 1:
                    rho = 1 / v1_pT(p, T);
                    break;
                case 2:
                    rho = 1 / v2_pT(p, T);
                    break;
                case 3:
                    rho = 1 / v3_ph(p, h3_pT(p, T));
                    break;
                case 4:
                    rho = 9999999999;  // error // need to check how to handle this
                    break;
                case 5:
                    rho = 1 / v5_pT(p, T);
                    break;
                default:
                    return 9999999999;  // error // need to check how to handle this
            }

            rhos = rho / 317.763;
            Ts = T / 647.226;
            ps = p / 22.115;

            //Check valid area
            if (T > 900 + 273.15 || (T > 600 + 273.15 && p > 300) || (T > 150 + 273.15 && p > 350) || p > 500)
            {
                return 9999999999;  // error // need to check how to handle this
            }
            my0 = Math.Pow(Ts, 0.5) / (1 + 0.978197 / Ts + 0.579829 / (Math.Pow(Ts, 2)) - 0.202354 / (Math.Pow(Ts, 3)));
            sum = 0;
            for (i = 0; i <= 5; i++)
            {
                sum = sum + h0[i] * Math.Pow(1 / Ts - 1, i) + h1[i] * Math.Pow(1 / Ts - 1, i) * Math.Pow(rhos - 1, 1) + h2[i] * Math.Pow(1 / Ts - 1, i) * Math.Pow(rhos - 1, 2) + h3[i] * Math.Pow(1 / Ts - 1, i) * Math.Pow(rhos - 1, 3) + h4[i] * Math.Pow(1 / Ts - 1, i) * Math.Pow(rhos - 1, 4) + h5[i] * Math.Pow(1 / Ts - 1, i) * Math.Pow(rhos - 1, 5) + h6[i] * Math.Pow(1 / Ts - 1, i) * Math.Pow(rhos - 1, 6);
            }
            my1 = Math.Exp(rhos * sum);
            return my0 * my1 * 0.000055071;
        }

        public double my_AllRegions_ph(double p, double h)
        {

            double rho = 0;
            double T = 0;
            double Ts = 0;
            double ps = 0;
            double my0 = 0;
            double sum = 0;
            double my1 = 0;
            double rhos = 0;
            double v4V = 0;
            double v4L = 0;
            double xs = 0;
            int i = 0;
            double[] h0 = { 0.5132047, 0.3205656, 0, 0, -0.7782567, 0.1885447 };
            double[] h1 = { 0.2151778, 0.7317883, 1.241044, 1.476783, 0, 0 };
            double[] h2 = { -0.2818107, -1.070786, -1.263184, 0, 0, 0 };
            double[] h3 = { 0.1778064, 0.460504, 0.2340379, -0.4924179, 0, 0 };
            double[] h4 = { -0.0417661, 0, 0, 0.1600435, 0, 0 };
            double[] h5 = { 0, -0.01578386, 0, 0, 0, 0 };
            double[] h6 = { 0, 0, 0, -0.003629481, 0, 0 };


            //Calcualte density.
            switch (region_ph(p, h))
            {
                case 1:
                    Ts = T1_ph(p, h);
                    T = Ts;
                    rho = 1 / v1_pT(p, Ts);
                    break;
                case 2:
                    Ts = T2_ph(p, h);
                    T = Ts;
                    rho = 1 / v2_pT(p, Ts);
                    break;
                case 3:
                    rho = 1 / v3_ph(p, h);
                    T = T3_ph(p, h);
                    break;
                case 4:
                    xs = x4_ph(p, h);
                    if (p < 16.529)
                    {
                        v4V = v2_pT(p, T4_p(p));
                        v4L = v1_pT(p, T4_p(p));
                    }
                    else
                    {
                        v4V = v3_ph(p, h4V_p(p));
                        v4L = v3_ph(p, h4L_p(p));
                    }
                    rho = 1 / (xs * v4V + (1 - xs) * v4L);
                    T = T4_p(p);
                    break;
                case 5:
                    Ts = T5_ph(p, h);
                    T = Ts;
                    rho = 1 / v5_pT(p, Ts);
                    break;
                default:
                    return 9999999999;  // error // need to check how to handle this
            }
            rhos = rho / 317.763;
            Ts = T / 647.226;
            ps = p / 22.115;
            //Check valid area
            if (T > 900 + 273.15 || (T > 600 + 273.15 && p > 300) || (T > 150 + 273.15 && p > 350) || p > 500)
            {
                return 9999999999;  // error // need to check how to handle this
            }
            my0 = Math.Pow(Ts, 0.5) / (1 + 0.978197 / Ts + 0.579829 / (Math.Pow(Ts, 2)) - 0.202354 / (Math.Pow(Ts, 3)));

            sum = 0;
            for (i = 0; i <= 5; i++)
            {
                sum = sum + h0[i] * Math.Pow(1 / Ts - 1, i) + h1[i] * Math.Pow(1 / Ts - 1, i) * Math.Pow(rhos - 1, 1) + h2[i] * Math.Pow(1 / Ts - 1, i) * Math.Pow(rhos - 1, 2) + h3[i] * Math.Pow(1 / Ts - 1, i) * Math.Pow(rhos - 1, 3) + h4[i] * Math.Pow(1 / Ts - 1, i) * Math.Pow(rhos - 1, 4) + h5[i] * Math.Pow(1 / Ts - 1, i) * Math.Pow(rhos - 1, 5) + h6[i] * Math.Pow(1 / Ts - 1, i) * Math.Pow(rhos - 1, 6);
            }
            my1 = Math.Exp(rhos * sum);
            return my0 * my1 * 0.000055071;
        }


        //***********************************************************************************************************
        //*5.2 Thermal Conductivity (IAPWS formulation 1985)

        public double tc_ptrho(double p, double T, double rho)
        {
            //Revised release on the IAPS Formulation 1985 for the Thermal Conductivity of ordinary water
            //IAPWS September 1998
            //Page 8
            //ver2.6 Start corrected bug
            double tc0 = 0;
            double tc1 = 0;
            double dT = 0;
            double Q = 0;
            double s = 0;
            double tc2 = 0;
         //   double tc = 0;
            if (T < 273.15)
            {
                return 9999999999;  // error // need to check how to handle this //Out of range of validity (para. B4)
            }
            else if (T < 500 + 273.15)
            {
                if (p > 100)
                {
                    return 9999999999;  // error // need to check how to handle this //Out of range of validity (para. B4)
                }
            }
            else if (T <= 650 + 273.15)
            {
                if (p > 70)
                {
                    return 9999999999;  // error // need to check how to handle this //Out of range of validity (para. B4)
                }
            }
            else if (T <= 800 + 273.15)
            {
                if (p > 40)
                {
                    return 9999999999;  // error // need to check how to handle this //Out of range of validity (para. B4)
                }
            }
            //ver2.6 End corrected bug

            T = T / 647.26;
            rho = rho / 317.7;
            tc0 = Math.Pow(T, 0.5) * (0.0102811 + 0.0299621 * T + 0.0156146 * Math.Pow(T, 2) - 0.00422464 * Math.Pow(T, 3));
            tc1 = -0.39707 + 0.400302 * rho + 1.06 * Math.Exp(-0.171587 * Math.Pow(rho + 2.39219, 2));
            dT = Math.Abs(T - 1) + 0.00308976;
            Q = 2 + 0.0822994 / Math.Pow(dT, 3 / 5.0);
            if (T >= 1)
            {
                s = 1 / dT;
            }
            else
            {
                s = 10.0932 / Math.Pow(dT, 3 / 5.0);
            }
            tc2 = (0.0701309 / Math.Pow(T, 10) + 0.011852) * Math.Pow(rho, 9 / 5.0) * Math.Exp(0.642857 * (1 - Math.Pow(rho, 14 / 5.0))) + 0.00169937 * s * Math.Pow(rho, Q) * Math.Exp((Q / (1 + Q)) * (1 - Math.Pow(rho, 1 + Q))) - 1.02 * Math.Exp(-4.11717 * Math.Pow(T, 3 / 2.0) - 6.17937 / Math.Pow(rho, 5));
            return tc0 + tc1 + tc2;
        }


        //***********************************************************************************************************
        //5.3 Surface Tension



        public double Surface_Tension_T(double T)
        {
            //IAPWS Release on Surface Tension of Ordinary Water Substance,
            //September 1994
            double tau = 0;
            double tc = 647.096;
            double b = 0.2358;
            double bb = -0.625;
            double my = 1.256;
            if (T < 0.01 || T > tc)
            {
                // return L"Out of valid region";
            }
            tau = 1 - T / tc;
            return b * Math.Pow(tau, my) * (1 + bb * tau);
        }

        //***********************************************************************************************************
        //*6 Units                                                                                      *
        //***********************************************************************************************************

        public double toSIunit_p(double Ins)
        {
            //Translate bar to MPa
            return Ins / 10;
        }

        public double fromSIunit_p(double Ins)
        {
            //Translate bar to MPa
            return Ins * 10;
        }

        public double toSIunit_T(double Ins)
        {
            //Translate degC to Kelvon
            return Ins + 273.15;
        }

        public double fromSIunit_T(double Ins)
        {
            //Translate Kelvin to degC
            return Ins - 273.15;
        }

        public double toSIunit_h(double Ins)
        {
            return Ins;
        }

        public double fromSIunit_h(double Ins)
        {
            return Ins;
        }

        public double toSIunit_v(double Ins)
        {
            return Ins;
        }

        public double fromSIunit_v(double Ins)
        {
            return Ins;
        }

        public double toSIunit_s(double Ins)
        {
            return Ins;
        }

        public double fromSIunit_s(double Ins)
        {
            return Ins;
        }

        public double toSIunit_u(double Ins)
        {
            return Ins;
        }

        public double fromSIunit_u(double Ins)
        {
            return Ins;
        }

        public double toSIunit_Cp(double Ins)
        {
            return Ins;
        }

        public double fromSIunit_Cp(double Ins)
        {
            return Ins;
        }

        public double toSIunit_Cv(double Ins)
        {
            return Ins;
        }

        public double fromSIunit_Cv(double Ins)
        {
            return Ins;
        }

        public double toSIunit_w(double Ins)
        {
            return Ins;
        }

        public double fromSIunit_w(double Ins)
        {
            return Ins;
        }

        public double toSIunit_tc(double Ins)
        {
            return Ins;
        }

        public double fromSIunit_tc(double Ins)
        {
            return Ins;
        }

        public double toSIunit_st(double Ins)
        {
            return Ins;
        }

        public double fromSIunit_st(double Ins)
        {
            return Ins;
        }

        public double toSIunit_x(double Ins)
        {
            return Ins;
        }

        public double fromSIunit_x(double Ins)
        {
            return Ins;
        }

        public double toSIunit_vx(double Ins)
        {
            return Ins;
        }

        public double fromSIunit_vx(double Ins)
        {
            return Ins;
        }

        public double toSIunit_my(double Ins)
        {
            return Ins;
        }

        public double fromSIunit_my(double Ins)
        {
            return Ins;
        }


    }

}


    

