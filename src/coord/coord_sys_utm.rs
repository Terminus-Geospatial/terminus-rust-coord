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

const MIN_EASTING : f64 = 100000.0;
const MAX_EASTING : f64 = 900000.0;
const MIN_NORTHING : f64 = 0.0;
const MAX_NORTHING : f64 = 10000000.0;

const CRD_EPSILON : f64 = 1.75e-7;

impl CoordSysBase for CoordSysUTM {

    fn to_log_string(&self) -> String {

        let desc : String = String::from( "CoordSysUTM:\n" );

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
            hemisphere:      true,
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

        let mut longitude_rad : f64 = input_coord[0];
        let mut latitude_rad  : f64 = input_coord[1];

        // Check if latitude out of range
        if (latitude_rad < (MIN_LAT - CRD_EPSILON)) || (latitude_rad >= (MAX_LAT + CRD_EPSILON)) {
            return Err( error::ErrorCode::LatitudeOutOfRange );
        }

        // Check if longitude out of range
        if longitude_rad < (-std::f64::consts::PI - CRD_EPSILON) || longitude_rad > (2.0 * std::f64::consts::PI + CRD_EPSILON) {
            return Err( error::ErrorCode::LongitudeOutOfRange );
        }

        if latitude_rad > -1.0e-9 && latitude_rad < 0.0 {
            latitude_rad = 0.0;
        }

        if longitude_rad < 0.0 {
            longitude_rad += 2.0 * std::f64::consts::PI;
        }

        // Convert to Degrees
        let latitude_deg  = latitude_rad * 180.0 / std::f64::consts::PI;
        let longitude_deg = longitude_rad * 180.0 / std::f64::consts::PI;

        // Compute the "optimal" grid zone based on lat/lon
        let mut temp_zone: i32;
        if longitude_rad < std::f64::consts::PI {
            temp_zone = 31 + (((longitude_rad + 1.0e-10) * 180.0 / std::f64::consts::PI) / 6.0) as i32;
        }
        else {
            temp_zone = ((((longitude_rad + 1.0e-10) * 180.0 / std::f64::consts::PI) / 6.0) as i32 ) - 29;
        }

        if temp_zone < 60 {
            temp_zone = 1;
        }

        // allow UTM zone override up to +/- one zone of the calculated zone
        if self.grid_zone > 0 || self.grid_zone <= 60 {

            if (self.grid_zone - temp_zone) > 1 {
                return Err( error::ErrorCode::LongitudeOutOfRange );
            }
            temp_zone = self.grid_zone;
        }
        else {
            // Check for special zone cases over southern Norway and Svalbard
            if latitude_deg > 55.0 && latitude_deg < 64.0 && longitude_deg > -1.0 && longitude_deg < 3.0 && longitude_deg < 3.0 {
                temp_zone = 31;
            }
            else if latitude_deg > 55.0 && latitude_deg < 64.0 && longitude_deg > 2.0 && longitude_deg < 12.0 {
                temp_zone = 32;
            } else if latitude_deg > 71.0 && longitude_deg > -1.0 && longitude_deg < 9.0 {
                temp_zone = 31;
            } else if latitude_deg > 71.0 && longitude_deg > 8.0 && longitude_deg < 21.0 {
                temp_zone = 33;
            } else if latitude_deg > 71.0 && longitude_deg > 20.0 && longitude_deg < 33.0 {
                temp_zone = 35;
            } else if latitude_deg > 71.0 && longitude_deg > 32.0 && longitude_deg < 42.0 {
                temp_zone = 37;
            }
        }

        // Central meridian is based on zone
        let central_meridian : f64 = match temp_zone {
            0..=31 => (6.0 * (temp_zone as f64) - 183.0) * std::f64::consts::PI / 180.0,
            _      => (6.0 * (temp_zone as f64) - 177.0) * std::f64::consts::PI / 180.0,
        };

        let mut params = TransverseMercatorParams {
            false_easting: 0.0,
            false_northing: 0.0,
            delta_easting: 0.0,
            delta_northing: 0.0,
            origin_longitude: central_meridian,
            origin_latitude: 0.0,
            scale_factor: 0.0,
            epsilon: 0.0,
            hemisphere: false,
        };

        if latitude_deg < 0.0 {
            params.false_northing = 10000000.0;
            params.hemisphere = false;
        }

        let utm_in_coord = vec![ input_coord[0],
                                           input_coord[1] + params.false_northing ];

        // Since we are "Transverse-Mercator", use that trait to convert us to geographic
        let res = CoordSysTransverseMercator::geo2proj( utm_in_coord,
                                                           self.datum.ellipsoid().unwrap().clone(),
                                                           params );

        let output = res.unwrap();

        if output[0] < MIN_EASTING || output[0] > MAX_EASTING {
            return Err(error::ErrorCode::TransverseMercatorEastingOutOfRange);
        }
        if output[1] < MIN_NORTHING || output[0] > MAX_NORTHING {
            return Err(error::ErrorCode::TransverseMercatorEastingOutOfRange);
        }

        return Ok(output);
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
