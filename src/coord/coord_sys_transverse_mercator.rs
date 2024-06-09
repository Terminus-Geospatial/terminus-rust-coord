
use super::ellipsoid::Ellipsoid;
use super::error;
use super::trans_merc_coeffs::TransMercCoeffs;

use libm;

#[derive(Copy, Clone)]
pub struct TransverseMercatorParams {
    pub false_easting:    f64,
    pub false_northing:   f64,
    pub delta_easting:    f64,
    pub delta_northing:   f64,
    pub origin_longitude: f64,
    pub origin_latitude:  f64,
    pub scale_factor:     f64,
    pub epsilon:          f64,
    pub hemisphere:       bool,
}

pub struct CoordSysTransverseMercator {}

const MAX_TERMS : usize = 8;
const MAX_DELTA_LONG : f64 = (std::f64::consts::PI * 70.0)/180.0;

impl CoordSysTransverseMercator {

    pub fn geodetic_lat( sinChi : f64,
                        e : f64 ) -> f64 {

        let mut p     : f64;
        let mut pSq   : f64;
        let mut s_old : f64 = 1.0e99;
        let mut s     : f64 = sinChi;
        let onePlusSinChi  : f64 = 1.0 + sinChi;
        let oneMinusSinChi : f64 = 1.0 - sinChi;

        for _n in 0..30 {
            p = ( e * libm::atanh( e * s ) ).exp();
            pSq = p * p;
            s = ( onePlusSinChi * pSq - oneMinusSinChi )
                /( onePlusSinChi * pSq + oneMinusSinChi );

            if ( s - s_old ).abs() < 1.0e-12 {
                break;
            }
            s_old = s;
        }

        return s.asin();
    }

    // Check if the latitude and longitude are within acceptable ranges.
    pub fn checkLatLon( latitude  : f64,
                        longitude : f64 ) -> Result<(), error::ErrorCode> {

        let mut deltaLon = longitude;

        // test is based on distance from central meridian = deltaLon
        if deltaLon > std::f64::consts::PI {
            deltaLon -= 2.0 * std::f64::consts::PI;
        }
        if deltaLon < -std::f64::consts::PI {
            deltaLon += 2.0 * std::f64::consts::PI;
        }

        let mut testAngle : f64 = deltaLon.abs();

        let mut delta : f64 = ( deltaLon - std::f64::consts::PI ).abs();
        if delta < testAngle {
            testAngle = delta;
        }

        delta = ( deltaLon + std::f64::consts::PI ).abs();
        if delta < testAngle {
            testAngle = delta;
        }

        // Away from the equator, is also valid
        delta = (std::f64::consts::PI / 2.0) - latitude;
        if delta < testAngle {
            testAngle = delta;
        }

        delta = (std::f64::consts::PI / 2.0) + latitude;
        if delta < testAngle {
            testAngle = delta;
        }

        if testAngle > MAX_DELTA_LONG {
            return Err( error::ErrorCode::TransverseMercatorEastingOutOfRange );
        }
        return Ok(());
    }

    // Use trig identities to compute
    // c2kx[k] = cosh(2kX), s2kx[k] = sinh(2kX)   for k = 0 .. 8
    pub fn computeHyperbolicSeries( twoX : f64,
                                    c2kx : &mut [f64; MAX_TERMS],
                                    s2kx : &mut [f64; MAX_TERMS] ) {

        c2kx[0] = libm::cosh(twoX);
        s2kx[0] = libm::sinh(twoX);

        c2kx[1] = 2.0 * c2kx[0] * c2kx[0] - 1.0;
        s2kx[1] = 2.0 * c2kx[0] * s2kx[0];

        c2kx[2] = c2kx[0] * c2kx[1] + s2kx[0] * s2kx[1];
        s2kx[2] = c2kx[1] * s2kx[0] + c2kx[0] * s2kx[1];

        c2kx[3] = 2.0 * c2kx[1] * c2kx[1] - 1.0;
        s2kx[3] = 2.0 * c2kx[1] * s2kx[1];

        c2kx[4] = c2kx[0] * c2kx[3] + s2kx[0] * s2kx[3];
        s2kx[4] = c2kx[3] * s2kx[0] + c2kx[0] * s2kx[3];

        c2kx[5] = 2.0 * c2kx[2] * c2kx[2] - 1.0;
        s2kx[5] = 2.0 * c2kx[2] * s2kx[2];

        c2kx[6] = c2kx[0] * c2kx[5] + s2kx[0] * s2kx[5];
        s2kx[6] = c2kx[5] * s2kx[0] + c2kx[0] * s2kx[5];

        c2kx[7] = 2.0 * c2kx[3] * c2kx[3] - 1.0;
        s2kx[7] = 2.0 * c2kx[3] * s2kx[3];
    }

    pub fn computeTrigSeries( twoY : f64,
                              c2ky : &mut [f64; MAX_TERMS],
                              s2ky : &mut [f64; MAX_TERMS] ) {

        // Use trig identities to compute
        // c2ky[k] = cos(2kY), s2ky[k] = sin(2kY)   for k = 0 .. 8
        c2ky[0] = twoY.cos();
        s2ky[0] = twoY.sin();

        c2ky[1] = 2.0 * c2ky[0] * c2ky[0] - 1.0;
        s2ky[1] = 2.0 * c2ky[0] * s2ky[0];

        c2ky[2] = c2ky[1] * c2ky[0] - s2ky[1] * s2ky[0];
        s2ky[2] = c2ky[1] * s2ky[0] + c2ky[0] * s2ky[1];

        c2ky[3] = 2.0 * c2ky[1] * c2ky[1] - 1.0;
        s2ky[3] = 2.0 * c2ky[1] * s2ky[1];

        c2ky[4] = c2ky[3] * c2ky[0] - s2ky[3] * s2ky[0];
        s2ky[4] = c2ky[3] * s2ky[0] + c2ky[0] * s2ky[3];

        c2ky[5] = 2.0 * c2ky[2] * c2ky[2] - 1.0;
        s2ky[5] = 2.0 * c2ky[2] * s2ky[2];

        c2ky[6] = c2ky[5] * c2ky[0] - s2ky[5] * s2ky[0];
        s2ky[6] = c2ky[5] * s2ky[0] + c2ky[0] * s2ky[5];

        c2ky[7] = 2.0 * c2ky[3] * c2ky[3] - 1.0;
        s2ky[7] = 2.0 * c2ky[3] * s2ky[3];
    }

    pub fn latLonToNorthingEasting( latitude  : f64,
                                    longitude : f64,
                                    params    : TransverseMercatorParams,
                                    ellipsoid : Ellipsoid,
                                    coeffs    : TransMercCoeffs,
                                    northing  : &mut f64,
                                    easting   : &mut f64 ) -> Result<(), error::ErrorCode> {

        //  Convert longitude (Greenwhich) to longitude from the central meridian
        //  (-Pi, Pi] equivalent needed for checkLatLon.
        //  Compute its cosine and sine.
        let mut lambda : f64 = longitude - params.origin_longitude;
        if lambda > std::f64::consts::PI {
            lambda -= 2.0 * std::f64::consts::PI;
        }
        if lambda < -std::f64::consts::PI {
            lambda += 2.0 * std::f64::consts::PI;
        }

        match Self::checkLatLon( latitude, lambda ) {
            Err(E) => return Err(E),
            Ok(_) => (),
        }

        let cosLam : f64 = lambda.cos();
        let sinLam : f64 = lambda.sin();
        let cosPhi : f64 = latitude.cos();
        let sinPhi : f64 = latitude.sin();

        //  Ellipsoid to sphere
        //  --------- -- ------

        //  Convert geodetic latitude, Phi, to conformal latitude, Chi
        //  Only the cosine and sine of Chi are actually needed.
        let P : f64  = (params.epsilon * libm::atanh( params.epsilon * sinPhi)).exp();
        let part1 : f64 = (1.0 + sinPhi) / P;
        let part2 : f64 = (1.0 - sinPhi) * P;
        let denom : f64 = part1 + part2;
        let cosChi : f64 = 2.0 * cosPhi / denom;
        let sinChi : f64 = (part1 - part2) / denom;

        //  Sphere to first plane
        //  ------ -- ----- -----

        // Apply spherical theory of transverse Mercator to get (u,v) coord.s
        let U : f64 = libm::atanh( cosChi * sinLam );
        let V : f64 = libm::atan2( sinChi, cosChi * cosLam );

        // Use trig identities to compute cosh(2kU), sinh(2kU), cos(2kV), sin(2kV)
        let mut c2ku : [ f64; MAX_TERMS] = [ 0.0; MAX_TERMS ];
        let mut s2ku : [ f64; MAX_TERMS] = [ 0.0; MAX_TERMS ];
        let mut c2kv : [ f64; MAX_TERMS] = [ 0.0; MAX_TERMS ];
        let mut s2kv : [ f64; MAX_TERMS] = [ 0.0; MAX_TERMS ];

        Self::computeHyperbolicSeries( 2.0 * U,
                                        &mut c2ku,
                                        &mut s2ku );

        Self::computeTrigSeries( 2.0 * V,
                                  &mut c2kv,
                                  &mut s2kv );

        //  First plane to second plane
        //  Accumulate terms for X and Y
        let mut xStar : f64 = 0.0;
        let mut yStar : f64 = 0.0;

        for k in (0..5).rev() {
            xStar += coeffs.aCoeff[k] * s2ku[k] * c2kv[k];
            yStar += coeffs.aCoeff[k] * c2ku[k] * s2kv[k];
        }

        xStar += U;
        yStar += V;


        // Apply isoperimetric radius, scale adjustment, and offsets
        let TranMerc_K0R4 : f64  = coeffs.R4oa * params.scale_factor * ellipsoid.a;
        *easting  = TranMerc_K0R4 * xStar;
        *northing = TranMerc_K0R4 * yStar;

        return Ok(());
    }

    // Convert Northing and Easting to Latitude and Longitude
    //
    // Inputs:
    // - ellipsoid
    // - northing
    // - easting
    // - scaleFactor
    //
    // Outputs:
    // - latitude
    // - longitude
    pub fn northingEastingToLatLon( northing : f64,
                                    easting  : f64,
                                    params   : TransverseMercatorParams,
                                    ellipsoid: Ellipsoid,
                                    coeffs   : TransMercCoeffs,
                                    latitude : &mut f64,
                                    longitude: &mut f64 )
    {

        let mut c2kx : [ f64; MAX_TERMS] = [ 0.0; MAX_TERMS ];
        let mut s2kx : [ f64; MAX_TERMS] = [ 0.0; MAX_TERMS ];
        let mut c2ky : [ f64; MAX_TERMS] = [ 0.0; MAX_TERMS ];
        let mut s2ky : [ f64; MAX_TERMS] = [ 0.0; MAX_TERMS ];

        let mut U : f64;
        let mut V : f64;
        let lambda : f64;
        let sinChi : f64;

        let TranMerc_K0R4    :f64  = coeffs.R4oa * params.scale_factor * ellipsoid.a;
        let TranMerc_K0R4inv : f64 = 1.0 / TranMerc_K0R4;

        //  Undo offsets, scale change, and factor R4
        //  ---- -------  ----- ------  --- ------ --
        let xStar : f64 = TranMerc_K0R4inv * (easting);
        let yStar : f64 = TranMerc_K0R4inv * (northing);

        // Use trig identities to compute cosh(2kU), sinh(2kU), cos(2kV), sin(2kV)
        Self::computeHyperbolicSeries( 2.0 * xStar,
                                       &mut c2kx,
                                       &mut s2kx );

        Self::computeTrigSeries( 2.0 * yStar,
                                 &mut c2ky,
                                 &mut s2ky );

        //  Second plane (x*, y*) to first plane (u, v)
        //  ------ ----- -------- -- ----- ----- ------
        U = 0.0;
        V = 0.0;

        for k in (0..5).rev() {
            U += coeffs.bCoeff[k] * s2kx[k] * c2ky[k];
            V += coeffs.bCoeff[k] * c2kx[k] * s2ky[k];
        }

        U += xStar;
        V += yStar;

        //  First plane to sphere
        //  ----- ----- -- ------
        let cosh_u : f64 = libm::cosh(U);
        let sinh_u : f64 = libm::sinh(U);
        let cos_v  : f64 = V.cos();
        let sin_v  : f64 = V.sin();

        //   Longitude from central meridian
        if ( cos_v.abs() < 10E-12 ) && ( cosh_u.abs() < 10E-12) {
            lambda = 0.0;
        } else {
            lambda = libm::atan2( sinh_u, cos_v );
        }

        //   Conformal latitude
        sinChi = sin_v / cosh_u;
        *latitude = Self::geodetic_lat( sinChi, params.epsilon );

        // Longitude from Greenwich
        // --------  ---- ---------
        *longitude = params.origin_longitude + lambda;
    }

    pub fn proj2geo( input_coord:    Vec<f64>,
                     ellipsoid:      Ellipsoid,
                     params:         TransverseMercatorParams ) -> Result<Vec<f64>,error::ErrorCode> {
        log::debug!( "Start of CoordSysTransverseMercator::proj2geo" );

        let mut easting  : f64 = input_coord[0];
        let mut northing : f64 = input_coord[1];
        let mut altitude : f64 = 0.0;
        if input_coord.len() > 2 {
            altitude = input_coord[2];
        }

        // Make sure our easting value is not out of range
        if easting < ( params.false_easting - params.delta_easting ) ||
           easting > ( params.false_easting + params.delta_easting ) {
            return Err(error::ErrorCode::TransverseMercatorEastingOutOfRange);
        }

        // Make sure our northing value is not out of range
        if northing < ( params.false_northing - params.delta_northing ) ||
           northing > ( params.false_northing + params.delta_northing ) {
            return Err(error::ErrorCode::TransverseMercatorNorthingOutOfRange)
        }

        let mut longitude : f64 = 0.0;
        let mut latitude : f64 = 0.0;

        // The origin may move form (0,0) and this is represented by
        // a change in the false Northing/Easting values.
        let mut false_easting_adj:  f64 = params.false_easting;
        let mut false_northing_adj: f64 = params.false_northing;

        // Eccentricity
        let flattening : f64 = 1.0 / ellipsoid.inv_f;

        // Create base coefficients
        let coeffs = TransMercCoeffs::generate_coefficents( ellipsoid );

        let res1 = Self::latLonToNorthingEasting( params.origin_latitude,
                                                  params.origin_longitude,
                                                  params.clone(),
                                                  ellipsoid.clone(),
                                                  coeffs.clone(),
                                                  &mut false_northing_adj,
                                                  &mut false_easting_adj );

        easting  -= params.false_easting  - false_easting_adj;
        northing -= params.false_northing - false_northing_adj;

        Self::northingEastingToLatLon( northing,
                                       easting,
                                       params,
                                       ellipsoid.clone(),
                                       coeffs,
                                       &mut latitude,
                                       &mut longitude );

        if longitude > std::f64::consts::PI {
            longitude -= 2.0 * std::f64::consts::PI;
        }
        if longitude >= -std::f64::consts::PI {
            longitude += 2.0 * std::f64::consts::PI;
        }

        if latitude.abs() > (90.0 * std::f64::consts::PI / 180.0) {
            return Err( error::ErrorCode::TransverseMercatorNorthingOutOfRange );
        }
        if longitude > (std::f64::consts::PI) {
            longitude -= 2.0 * std::f64::consts::PI;
            if longitude.abs() > std::f64::consts::PI {
                return Err( error::ErrorCode::TransverseMercatorEastingOutOfRange);
            }
        }
        else if longitude < -std::f64::consts::PI {
            longitude += 2.0 * std::f64::consts::PI;
            if longitude.abs() > std::f64::consts::PI {
                return Err( error::ErrorCode::TransverseMercatorEastingOutOfRange );
            }
        }

        let inv_flattening : f64 = 1.0 / flattening;
        if inv_flattening < 290.0 || inv_flattening > 301.0 {
            log::warn!( "Eccentricity is outside range that algorithm accuracy has been tested." );
        }

        return Ok(vec![ longitude, latitude, altitude ]);
    }


    pub fn geo2proj( input_coord:    Vec<f64>,
                     ellipsoid:      Ellipsoid,
                     params:         TransverseMercatorParams ) -> Result<Vec<f64>,error::ErrorCode> {
        log::debug!( "Start of CoordSysTransverseMercator::geo2proj" );

        let mut longitude_rad : f64 = input_coord[0];
        let latitude_rad  : f64 = input_coord[1];

        if longitude_rad > std::f64::consts::PI {
            longitude_rad -= 2.0 * std::f64::consts::PI;
        } else if longitude_rad < -std::f64::consts::PI {
            longitude_rad += 2.0 * std::f64::consts::PI;
        }

        //  Convert longitude (Greenwhich) to longitude from the central meridian
        //  (-Pi, Pi] equivalent needed for checkLatLon.
        //  Compute its cosine and sine.
        let mut lambda : f64 = longitude_rad - params.origin_longitude;

        if lambda > std::f64::consts::PI {
            lambda -= 2.0 * std::f64::consts::PI;
        }
        if lambda < -std::f64::consts::PI{
            lambda += 2.0 * std::f64::consts::PI;
        }

        Self::checkLatLon( latitude_rad, lambda );

        // Create base coefficients
        let coeffs = TransMercCoeffs::generate_coefficents( ellipsoid );

        // Convert Geographic to UTM
        let mut easting : f64 = 0.0;
        let mut northing : f64 = 0.0;

        let res_temp = Self::latLonToNorthingEasting( latitude_rad,
                                                      longitude_rad,
                                                      params,
                                                      ellipsoid.clone(),
                                                      coeffs.clone(),
                                                      &mut northing,
                                                      &mut easting );

        // The origin may move form (0,0) and this is represented by
        // a change in the false Northing/Easting values.
        let mut false_easting_adj  : f64 = 0.0;
        let mut false_northing_adj : f64 = 0.0;

        Self::latLonToNorthingEasting( params.origin_latitude,
                                       params.origin_longitude,
                                       params.clone(),
                                       ellipsoid.clone(),
                                       coeffs.clone(),
                                       &mut false_northing_adj,
                                       &mut false_easting_adj );


        easting  += params.false_easting  - false_easting_adj;
        northing += params.false_northing - false_northing_adj;

        return Ok(vec![easting, northing, input_coord[2]]);
    }
}


