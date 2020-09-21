#pragma once
#include<iostream>

namespace mswhgeo
{
    
struct GeoEllipsoid
{
    double a;
    double f_inv;
    
    double f;
    double b;
    double e;
    double e_2;
    double hmin;
    
    GeoEllipsoid(double a_, double f_inv_);
    GeoEllipsoid(const GeoEllipsoid& gE);
};   

class Geo
{
public:
    Geo(GeoEllipsoid gE);
    void Wgs2Ecef(double lat, double lon, double alt, double& x, double& y, double& z);
    int Ecef2Wgs(double x, double y, double z, double& lat, double& lon, double& alt);
        
    void Wgs2Enu(/*ref point*/double lat0, double lon0, double alt0, /*wgs point*/double lat, double lon, double alt, /*out enu point*/ double& e, double& n, double& u);
    int Enu2Wgs(/*ref point*/double lat0, double lon0, double alt0, /*enu point*/double e, double n, double u, /*out wgs point*/ double& lat, double& lon, double& alt);
    
    void Ecef2Enu(/*ref point*/double lat0, double lon0, double alt0, /*ecef point*/double x, double y, double z,  /*out enu point*/double& e, double& n, double& u);
    void Enu2Ecef(/*ref point*/double lat0, double lon0, double alt0, /*enu point*/double e, double n, double u,  /*out ecef point*/double& x, double& y, double& z);
    
private:
    GeoEllipsoid _gE;
};

}

namespace mswhgeo_constants{
    static const struct mswhgeo::GeoEllipsoid WGS84 = {6378137.0, 298.257223563};
    static const double PI = 3.1415926535898;
}
