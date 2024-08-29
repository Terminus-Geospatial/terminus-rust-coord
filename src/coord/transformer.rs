
use super::coord_sys_base::CoordSysBase;
use super::error;
use std::vec::Vec;

pub struct Transformer {
    pub cs_in:  Box<dyn CoordSysBase>,
    pub cs_out: Box<dyn CoordSysBase>,
}

impl Transformer {

    pub fn create( cs_in: Box<dyn CoordSysBase>, cs_out: Box<dyn CoordSysBase> ) -> Box<Transformer> {
        return Box::new( Transformer {
            cs_in:  cs_in,
            cs_out: cs_out
        } );
    }


    pub fn apply(&self, input_coord: Vec<f64> ) -> Result<Vec<f64>,error::ErrorCode> {

        // Step 1:  Convert input CS to geographic
        let geo_res = match self.cs_in.is_projected() {
            true => {
                // Handle all errors if the call actually fails
                match self.cs_in.to_geographic( input_coord ) {
                    Ok(t) => t,
                    Err(e) => return Err(e),
                }
            },
            false => return Ok(input_coord),
        };

        // Step 2:  Perform Geographic-to-Geographic (TODO)
        let geo_final = geo_res;

        // Step 3: Convert to final projected
        let proj_res = match self.cs_out.is_projected() {
            true => {
                match self.cs_out.from_geographic( geo_final ) {
                    Ok(t) => t,
                    Err( e ) => return Err(e),
                }
            },
            false => geo_final,
        };

        Ok(proj_res)
    }
}
