#include<iostream>
#include<mshw_geolib.h>
#include<cmath>

double subtractAngles(double x, double y)
{
    double pi = mswhgeo_constants::PI;
    return fmin(fabs(x - y), fabs(fabs(x - y) - 2 * pi));
}
 
int main(int argc, char *argv[])
{
   double meanAngleError = 0.0;
   double meanAltError = 0.0;
   double meanEnuError = 0.0;
   int N = 0;
   int K = 0;
   double pi = mswhgeo_constants::PI;
   mswhgeo::Geo geo(mswhgeo_constants::WGS84);
   for (double lat = -85 * pi / 180.0; lat < 85 * pi / 180.0; lat += 10.0 * pi / 180.0)
   {
   for (double lon = -180 * pi / 180.0; lon < 180 * pi / 180.0; lon += 10.0 * pi / 180.0)
   {
   for (double alt = -5000.0; alt < 5000.0; alt += 1000.0)
   {
        N++;
       
        double x = 0;
        double y = 0;
        double z = 0;
       
        double lat_tr = 0;
        double lon_tr = 0;
        double alt_tr = 0;
        geo.Wgs2Ecef(lat, lon, alt, x, y, z);
        geo.Ecef2Wgs(x, y, z, lat_tr, lon_tr, alt_tr);
        
        meanAngleError += subtractAngles(lat, lat_tr) + subtractAngles(lon, lon_tr);
        meanAltError += fabs(alt_tr - alt);
        
        for (double dlat = -0.1 * pi / 180.0; dlat < 0.1 * pi / 180.0; dlat += 0.01 * pi / 180.0)
        {
        for (double dlon = -0.1 * pi / 180.0; dlon < 0.1 * pi / 180.0; dlon += 0.01 * pi / 180.0)
        {
        for (double dalt = -500.0; dalt < 500.0; dalt += 500.0)
        {
            
        K++;
        double e = 0;
        double n = 0;
        double u = 0;
        double e_tr = 0;
        double n_tr = 0;
        double u_tr = 0;
        geo.Wgs2Enu(lat, lon, alt, lat + dlat, lon + dlon, alt + dalt, e, n, u);
        geo.Enu2Wgs(lat, lon, alt, e, n, u, lat_tr, lon_tr, alt_tr);
        geo.Wgs2Enu(lat, lon, alt, lat_tr, lon_tr, alt_tr, e_tr, n_tr, u_tr);
        
        meanEnuError += fabs(e - e_tr) + fabs(n - n_tr) +fabs(u - u_tr);

        }
        }
        }
   }
   }
   }
   
    printf("meanAngleError: %.40f\n", meanAngleError/N/2.0);
    printf("meanAltError: %.40f\n", meanAltError/N);
    printf("meanEnuError: %.40f\n", meanEnuError/K/3.0);
    return 0;
}
