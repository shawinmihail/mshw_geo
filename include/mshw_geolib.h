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
    
    GeoEllipsoid(double a_, double f_inv_);
    GeoEllipsoid(const GeoEllipsoid& gE);
};   

class Geo
{
public:
    Geo(GeoEllipsoid gE);
    int Wgs2Ecef(double lat, double lon, double alt, double& x, double& y, double& z);
    int Ecef2Wgs(double x, double y, double z, double& lat, double& lon, double& alt);
        
    int Wgs2Enu(double lat, double lon, double alt, double& x, double& y, double& z);
    int Enu2Wgs(double lat, double lon, double alt, double& x, double& y, double& z);
    
    void rotateEcef2Enu(/*ref point*/double lat0, double lon0, double alt0, /*ecef point*/double x, double y, double z,  /*out enu point*/double& e, double& n, double& u);
    void rotateEnu2Ecef(/*ref point*/double lat0, double lon0, double alt0, /*enu point*/double e, double n, double u,  /*out ecef point*/double& x, double& y, double& z);
    
private:
    GeoEllipsoid _gE;
};

}

namespace mswhgeo_constants{
    static const struct mswhgeo::GeoEllipsoid WGS84 = {6378137.0, 298.257223563};
    static const double PI = 3.1415926535898;
}
