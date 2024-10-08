/**************************** INTELLECTUAL PROPERTY RIGHTS ****************************/
/*                                                                                    */
/*                           Copyright (c) 2024 Terminus LLC                          */
/*                                                                                    */
/*                                All Rights Reserved.                                */
/*                                                                                    */
/*          Use of this source code is governed by LICENSE in the repo root.          */
/*                                                                                    */
/***************************# INTELLECTUAL PROPERTY RIGHTS ****************************/

// Do not delete this EVER!
//
// From GeoTrans 3.9, ellips.dat
//
//Airy 1830                      AA 6377563.396000000 6356256.909237285 299.3249646000000
//Modified Airy                  AM 6377340.189000000 6356034.447938534 299.3249646000000
//Australian National            AN 6378160.000000000 6356774.719195306 298.2500000000000
//Bessel 1841 (Namibia)          BN 6377483.865000000 6356165.382966325 299.1528128000000
//Bessel 1841 (Ethiopia, etc.)   BR 6377397.155000000 6356078.962818188 299.1528128000000
//Clarke 1866                    CC 6378206.400000000 6356583.800000000 294.9786982139058
//Clarke 1880                    CD 6378249.145000000 6356514.869549776 293.4650000000000
//Clarke 1880 (IGN)              CG 6378249.200000000 6356514.999963442 293.4660208000000
//Everest (India 1830)           EA 6377276.345000000 6356075.413140240 300.8017000000000
//Everest (E. Malaysia, Brunei)  EB 6377298.556000000 6356097.550300897 300.8017000000000
//Everest (India 1956)           EC 6377301.243000000 6356100.228368101 300.8017000000000
//Everest (W. Malaysia 1969)     ED 6377295.664000000 6356094.667915204 300.8017000000000
//Everest (W. Mal. & Sing. 1948) EE 6377304.063000000 6356103.038993154 300.8017000000000
//Everest (Pakistan)             EF 6377309.613000000 6356108.570542461 300.8017000000000
//Mod. Fischer 1960 (S. Asia)    FA 6378155.000000000 6356773.320482736 298.3000000000000
//Helmert 1906                   HE 6378200.000000000 6356818.169627891 298.3000000000000
//Hough 1960                     HO 6378270.000000000 6356794.343434343 297.0000000000000
//Indonesian 1974                ID 6378160.000000000 6356774.504085540 298.2470000000000
//International 1924             IN 6378388.000000000 6356911.946127946 297.0000000000000
//Krassovsky 1940                KA 6378245.000000000 6356863.018773047 298.3000000000000
//GRS 80                         RF 6378137.000000000 6356752.314140356 298.2572221010000
//South American 1969            SA 6378160.000000000 6356774.719195306 298.2500000000000
//WGS 72                         WD 6378135.000000000 6356750.520016093 298.2600000000000
//WGS 84                         WE 6378137.000000000 6356752.314245179 298.2572235630000
//War Office 1924                WO 6378300.580000000 6356752.267229729 296.0000000000000

#[derive(Copy, Clone)]
pub enum EllipsoidType {
    CUSTOM,
    WGS84,
}

#[derive(Copy, Clone)]
pub struct Ellipsoid {
    pub a: f64,
    pub b: f64,
    pub inv_f: f64,
    pub code: EllipsoidType,
}

impl Ellipsoid {

    pub fn new( a: f64, b: f64, inv_f: f64, code: EllipsoidType ) -> Ellipsoid {
        let d = Ellipsoid { 
            a:      a,
            b:      b,
            inv_f: inv_f,
            code:  code
        };
        return d;
    }

    pub fn wgs84() -> Ellipsoid {
        Ellipsoid::new( 6378137.0,
                        6356752.314245,
                        298.257223563,
                        EllipsoidType::WGS84 )
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_create_wgs84() {
        let wgs84 = Ellipsoid::wgs84();
        assert_eq!( ( wgs84.a - 6378137.0 ).abs() < 0.00001, true );
    }
}