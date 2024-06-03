use crate::coord::coord_sys_def::CoordSysDefinition;
use crate::coord::{datum, error};
use super::coord_sys_base::CoordSysBase;
use super::coord_sys_transverse_mercator::{CoordSysTransverseMercator, TransverseMercatorParams};

pub struct CoordSysUTM {
    datum: Box<dyn datum::Datum>,
    grid_zone: i32,
    is_northern_hemi: bool,
}

const MIN_LAT : f64 = (-80.5 * std::f64::consts::PI) / 180.0; // -80.5 degrees in radians
const MAX_LAT : f64 = (84.5 * std::f64::consts::PI) / 180.0;  //  84.5 degrees in radians

const CRD_EPSILON : f64 = 1.75e-7;

impl CoordSysBase for CoordSysUTM {

    fn to_log_string(&self) -> String {

        let mut desc : String = String::from( "CoordSysUTM:\n" );

        return desc;
    }

    fn to_geographic(&self, input_coord: Vec<f64>) -> Result<Vec<f64>,error::ErrorCode> {
        log::debug!( "Start of CoordSysUTM::to_geographic" );

        // Central meridian is based on zone
        let central_meridian : f64 = match self.grid_zone {
            0..=31 => (6.0 * (self.grid_zone as f64) - 183.0) * std::f64::consts::PI / 180.0,
            _      => (6.0 * (self.grid_zone as f64) - 177.0) * std::f64::consts::PI / 180.0,
        };

        let flattening = 1.0 / self.datum.ellipsoid().unwrap().inv_f;

        // Constants for method
        let mut params = TransverseMercatorParams {
            false_easting:    500000.0,
            false_northing:   0.0,
            delta_easting:    10000000.0,
            delta_northing:   20000000.0,
            origin_longitude: central_meridian,
            origin_latitude:  0.0,
            scale_factor:     0.9996,
            epsilon: ( 2.0 * flattening - flattening * flattening ).sqrt(),
        };

        // False Northing depends on hemisphere
        if !self.is_northern_hemi {
            params.false_northing = 10000000.0;
        }
        
        // Since we are "Transverse-Mercator", use that trait to convert us to geographic
        return CoordSysTransverseMercator::proj2geo( input_coord,
                                                     self.datum.ellipsoid().unwrap().clone(),
                                                     params );
    }

    fn from_geographic(&self, input_coord: Vec<f64>) -> Result<Vec<f64>,error::ErrorCode> {
        log::debug!( "Start of CoordSysUTM::from_geographic" );

        let longitude : f64 = input_coord[0];
        let latitude  : f64 = input_coord[1];
        let mut altitude  : f64 = 0.0;
        if input_coord.len() > 2 {
            altitude = input_coord[2];
        }


        long Lat_Degrees;
        long Long_Degrees;
        long temp_zone;
        char hemisphere;
        double False_Northing = 0;


        // Check if latitude out of range
        if (latitude < (MIN_LAT - CRD_EPSILON)) || (latitude >= (MAX_LAT + CRD_EPSILON)) {
            return Err( error::ErrorCode::LatitudeOutOfRange );
        }

        // Check if longitude out of range
        if longitude < (-std::f64::consts::PI - CRD_EPSILON) || longitude > (2*std::f64::consts::PI + CRD_EPSILON) {
            return Err( error::ErrorCode::LongitudeOutOfRange );
        }

        if((latitude > -1.0e-9) && (latitude < 0))
        latitude = 0.0;

        if (longitude < 0)
        longitude += (2*PI);

        Lat_Degrees = (long)(latitude * 180.0 / PI);
        Long_Degrees = (long)(longitude * 180.0 / PI);

        if (longitude < PI)
        temp_zone = (long)(31 + (((longitude+1.0e-10) * 180.0 / PI) / 6.0));
        else
        temp_zone = (long)((((longitude+1.0e-10) * 180.0 / PI) / 6.0) - 29);

        if (temp_zone > 60)
        temp_zone = 1;

        /* allow UTM zone override up to +/- one zone of the calculated zone */
        if( utmZoneOverride )
        {
            if ((temp_zone == 1) && (utmZoneOverride == 60))
            temp_zone = utmZoneOverride;
            else if ((temp_zone == 60) && (utmZoneOverride == 1))
            temp_zone = utmZoneOverride;
            else if (((temp_zone-1) <= utmZoneOverride) &&
            (utmZoneOverride <= (temp_zone+1)))
            temp_zone = utmZoneOverride;
            else
            throw CoordinateConversionException( ErrorMessages::zoneOverride );
        }
        else if( UTM_Override )
        {
            if ((temp_zone == 1) && (UTM_Override == 60))
            temp_zone = UTM_Override;
            else if ((temp_zone == 60) && (UTM_Override == 1))
            temp_zone = UTM_Override;
            else if (((temp_zone-1) <= UTM_Override) &&
            (UTM_Override <= (temp_zone+1)))
            temp_zone = UTM_Override;
            else
            throw CoordinateConversionException( ErrorMessages::zoneOverride );
        }
        else /* not UTM zone override */
        {
            /* check for special zone cases over southern Norway and Svalbard */
            if ((Lat_Degrees > 55) && (Lat_Degrees < 64) && (Long_Degrees > -1)
                && (Long_Degrees < 3))
            temp_zone = 31;
            if ((Lat_Degrees > 55) && (Lat_Degrees < 64) && (Long_Degrees > 2)
                && (Long_Degrees < 12))
            temp_zone = 32;
            if ((Lat_Degrees > 71) && (Long_Degrees > -1) && (Long_Degrees < 9))
            temp_zone = 31;
            if ((Lat_Degrees > 71) && (Long_Degrees > 8) && (Long_Degrees < 21))
            temp_zone = 33;
            if ((Lat_Degrees > 71) && (Long_Degrees > 20) && (Long_Degrees < 33))
            temp_zone = 35;
            if ((Lat_Degrees > 71) && (Long_Degrees > 32) && (Long_Degrees < 42))
            temp_zone = 37;
        }

        TransverseMercator *transverseMercator = transverseMercatorMap[temp_zone];

        if (latitude < 0)
        {
            False_Northing = 10000000;
            hemisphere = 'S';
        }
        else
        hemisphere = 'N';

        GeodeticCoordinates tempGeodeticCoordinates(
            CoordinateType::geodetic, longitude, latitude );
        MapProjectionCoordinates* transverseMercatorCoordinates =
            transverseMercator->convertFromGeodetic( &tempGeodeticCoordinates );
        double easting = transverseMercatorCoordinates->easting();
        double northing = transverseMercatorCoordinates->northing() + False_Northing;

        if ((easting < MIN_EASTING) || (easting > MAX_EASTING))
        {
            delete transverseMercatorCoordinates;
            throw CoordinateConversionException( ErrorMessages::easting );
        }

        if ((northing < MIN_NORTHING) || (northing > MAX_NORTHING))
        {
            delete transverseMercatorCoordinates;
            throw CoordinateConversionException( ErrorMessages::northing );
        }

        delete transverseMercatorCoordinates;

        return new UTMCoordinates(
            CoordinateType::universalTransverseMercator,
            temp_zone, hemisphere, easting, northing );

        Ok(input_coord.to_vec())
    }

    fn is_geographic(&self) -> bool {
        false
    }

    fn is_projected(&self) -> bool {
        true
    }
}


impl CoordSysUTM {

    pub fn create( cs_def : CoordSysDefinition ) -> Option<Box<dyn CoordSysBase>> {

        // If the EPSG code is specified, use this to build our cs
        match cs_def.epsg_code {
            Some(x) => {

                // WGS84, standard UTM
                match x {
                    32600..=32700 => {
                        return Some( Box::new( CoordSysUTM {
                            datum: Box::new( datum::wgs84::WGS84::as_epsg_4326() ),
                            grid_zone: x - 32600,
                            is_northern_hemi: true,
                        }));
                    },
                    _ => return None,
                }
            },
            None => (), // do nothing
        }

        return None;
    }
}
