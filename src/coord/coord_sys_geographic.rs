/**************************** INTELLECTUAL PROPERTY RIGHTS ****************************/
/*                                                                                    */
/*                           Copyright (c) 2024 Terminus LLC                          */
/*                                                                                    */
/*                                All Rights Reserved.                                */
/*                                                                                    */
/*          Use of this source code is governed by LICENSE in the repo root.          */
/*                                                                                    */
/***************************# INTELLECTUAL PROPERTY RIGHTS ****************************/

use super::coord_sys_def::CoordSysDefinition;
use super::coord_sys_base::CoordSysBase;
use super::datum;
use super::error;

pub struct CoordSysGeographic {
    datum: Box<dyn datum::Datum>,
}

impl CoordSysBase for CoordSysGeographic {

    fn to_log_string(&self) -> String {

        let desc : String = String::from( "CoordSysGeographic:\n" );

        return desc;
    }

    fn to_geographic(&self, input_coord: Vec<f64>) -> Result<Vec<f64>,error::ErrorCode> {

        log::debug!( "Start of CoordSysGeographic::to_geographic" );
        Ok(input_coord.to_vec())
    }

    fn from_geographic(&self, input_coord: Vec<f64>) -> Result<Vec<f64>,error::ErrorCode> {
        log::debug!( "Start of CoordSysGeographic::from_geographic" );
        Ok(input_coord.to_vec())
    }

    fn is_geographic(&self) -> bool {
        true
    }

    fn is_projected(&self) -> bool {
        false
    }
}

impl CoordSysGeographic {

    pub fn create( cs_def : CoordSysDefinition ) -> Option<Box<dyn CoordSysBase>> {

        // If the EPSG is 4326, this is really straightforward
        match cs_def.epsg_code {
            Some(4326) => {
                return Some( Box::new( CoordSysGeographic {
                    datum: Box::new( datum::wgs84::WGS84::as_epsg_4326() ),
                } ) );
            },
            Some(_) => return None,
            None => (),
        }
        return None;
    }
}